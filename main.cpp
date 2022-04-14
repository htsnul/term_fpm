#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <memory>
#include <vector>

#include <signal.h>
#include <stdlib.h>
#include <sys/select.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>

class Terminal {
public:
  inline static bool exitRequested{};
  static void initialize() {
    {
      struct sigaction sa {};
      sa.sa_handler = [](int) { exitRequested = true; };
      sigaction(SIGINT, &sa, nullptr);
    }
    // https://www.gnu.org/software/libc/manual/html_node/Noncanon-Example.html
    {
      tcgetattr(STDIN_FILENO, &originalTermios_);
      termios t{originalTermios_};
      t.c_lflag &= ~(ICANON | ECHO);
      t.c_cc[VMIN] = 1;
      t.c_cc[VTIME] = 0;
      tcsetattr(STDOUT_FILENO, TCSAFLUSH, &t);
    }
    hideCursor();
    eraseEntireScreen();
  }

  static void finalize() {
    showCursor();
    tcsetattr(STDOUT_FILENO, TCSANOW, &originalTermios_);
  }

  static void eraseEntireScreen() { std::printf("\x1b[2J"); }

  static void setCursorPositionToTopLeft() { std::printf("\x1b[1;1H"); }

  static bool canReadStdin() {
    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(STDIN_FILENO, &fds);
    timeval tv{};
    return select(1, &fds, nullptr, nullptr, &tv) > 0;
  }

private:
  static void hideCursor() { std::printf("\x1b[?25l"); }

  static void showCursor() { std::printf("\x1b[?25h"); }

  inline static termios originalTermios_{};
};

int64_t getCurrentTimeMs() {
  struct timespec ts {};
  clock_gettime(CLOCK_REALTIME, &ts);
  return 1000 * int64_t{ts.tv_sec} + int64_t{ts.tv_nsec} / (1000 * 1000);
}

void sleepTillMs(int64_t timeMs) {
  const int64_t currentTimeMs = getCurrentTimeMs();
  if (timeMs <= currentTimeMs)
    return;
  usleep(1000 * (timeMs - currentTimeMs));
}

struct Pixel {
  uint8_t r, g, b;
};

class Canvas {
public:
  inline static const int WIDTH{80}, HEIGHT{40};
  void clear();
  Pixel &pixel(int x, int y) { return pixels[y][x]; }
  void printToTerminal();

private:
  Pixel pixels[HEIGHT][WIDTH];
};

void Canvas::clear() {
  for (int y = 0; y < HEIGHT; ++y) {
    for (int x = 0; x < WIDTH; ++x) {
      pixels[y][x] = {0, 0, 0};
    }
  }
}

void Canvas::printToTerminal() {
  Terminal::setCursorPositionToTopLeft();
  for (int hy = 0; hy < HEIGHT / 2; ++hy) {
    for (int x = 0; x < WIDTH; ++x) {
      const auto &pu = pixels[2 * hy + 0][x];
      const auto &pd = pixels[2 * hy + 1][x];
      const int dou = (int[]){+5, -5}[x % 2];
      const int dod = (int[]){-15, +15}[x % 2];
      const Pixel dpu{
          static_cast<uint8_t>(std::clamp(pu.r + dou, 0, 255)),
          static_cast<uint8_t>(std::clamp(pu.g + dou, 0, 255)),
          static_cast<uint8_t>(std::clamp(pu.b + dou, 0, 255)),
      };
      const Pixel dpd{
          static_cast<uint8_t>(std::clamp(pd.r + dod, 0, 255)),
          static_cast<uint8_t>(std::clamp(pd.g + dod, 0, 255)),
          static_cast<uint8_t>(std::clamp(pd.b + dod, 0, 255)),
      };
      std::printf("\x1b[38;5;%dm\x1b[48;5;%dm\u2580",
                  16 + 36 * (dpu.r * 6 / 256) + 6 * (dpu.g * 6 / 256) +
                      (dpu.b * 6 / 256),
                  16 + 36 * (dpd.r * 6 / 256) + 6 * (dpd.g * 6 / 256) +
                      (dpd.b * 6 / 256));
    }
    std::printf("\n");
  }
}

class InputManager {
public:
  static bool key(char c) { return keys_[c]; }

  static void update() {
    char data[256]{};
    ssize_t readSize{};
    if (Terminal::canReadStdin()) {
      readSize = read(STDIN_FILENO, data, std::size(data));
    }
    for (int i = 0; i < 256; ++i) {
      keys_[i] = false;
    }
    for (int i = 0; i < readSize; ++i) {
      const char c = data[i];
      keys_[c] = true;
    }
  }

  inline static bool keys_[256]{};
};

struct Vec3 {
  float x, y, z;
  float length() const { return std::sqrt(x * x + y * y + z * z); }
  Vec3 add(const Vec3 &v) const { return {x + v.x, y + v.y, z + v.z}; }
  Vec3 sub(const Vec3 &v) const { return add({-v.x, -v.y, -v.z}); }
  Vec3 mul(float s) const { return {s * x, s * y, s * z}; }
  float dot(const Vec3 &v) const { return {x * v.x + y * v.y + z * v.z}; }
  Vec3 mod(float s) const {
    return {
        std::fmod(x, s) + (x < 0 ? s : 0),
        std::fmod(y, s) + (y < 0 ? s : 0),
        std::fmod(z, s) + (z < 0 ? s : 0),
    };
  }
  Vec3 normalized() {
    const float l = length();
    return l == 0.0f ? Vec3{} : mul(1.0f / l);
  }
  Vec3 abs() const {
    return {
        std::abs(x),
        std::abs(y),
        std::abs(z),
    };
  }
  Vec3 max(float lo) const {
    return {
        std::max(x, lo),
        std::max(y, lo),
        std::max(z, lo),
    };
  }
  Vec3 clamped(float lo, float hi) const {
    return {
        std::clamp(x, lo, hi),
        std::clamp(y, lo, hi),
        std::clamp(z, lo, hi),
    };
  }
  Vec3 rotateY(float angle) const {
    return {
        std::cos(angle) * x - std::sin(angle) * z,
        y,
        std::sin(angle) * x + std::cos(angle) * z,
    };
  }
};

class Maze {
public:
  inline static const int WIDTH{9}, DEPTH{9};
  void initialize() {
    for (int z = 0; z < DEPTH; ++z) {
      for (int x = 0; x < WIDTH; ++x) {
        blocks[z][x] = true;
      }
    }
    recurse(1, 1);
  }
  void recurse(int z, int x) {
    blocks[z][x] = false;
    while (true) {
      auto neighbors = getUnvisitedNeighbors(z, x);
      if (neighbors.empty()) return;
      auto& neighbor = neighbors.at(std::rand() % neighbors.size());
      const int nz = neighbor.first;
      const int nx = neighbor.second;
      blocks[(z + nz) / 2][(x + nx) / 2] = false;
      recurse(nz, nx);
    }
  }
  std::vector<std::pair<int, int>> getUnvisitedNeighbors(int z, int x) {
    std::vector<std::pair<int, int>> neighbors{};
    if (z - 2 >= 0 && blocks[z - 2][x])
      neighbors.push_back(std::pair<int, int>(z - 2, x));
    if (x - 2 >= 0 && blocks[z][x - 2])
      neighbors.push_back(std::pair<int, int>(z, x - 2));
    if (z + 2 < DEPTH && blocks[z + 2][x])
      neighbors.push_back(std::pair<int, int>(z + 2, x));
    if (x + 2 < WIDTH && blocks[z][x + 2])
      neighbors.push_back(std::pair<int, int>(z, x + 2));
    return neighbors;
  }
  float getDistance(const Vec3 &pos) const {
    const int sx = std::max(0, static_cast<int>(std::round(pos.x)) - 1);
    const int ex = std::min(WIDTH, sx + 2);
    const int sz = std::max(0, static_cast<int>(std::round(pos.z)) - 1);
    const int ez = std::min(DEPTH, sz + 2);
    float m = std::min({pos.y - (-0.5f), 0.5f - pos.y});
    for (int z = sz; z < ez; ++z) {
      for (int x = sx; x < ex; ++x) {
        if (!blocks[z][x])
          continue;
        const Vec3 bp{static_cast<float>(x) + 0.5f, 0.0f,
                      static_cast<float>(z) + 0.5f};
        m = std::min(
            m,
            pos.sub(bp).abs().sub(Vec3{0.5f, 0.5f, 0.5f}).max(0.0f).length());
      }
    }
    return m;
  }

  Vec3 getNormal(const Vec3 &pos) const {
    const float ep = 0.0001f;
    const Vec3 v0{pos.x - ep, pos.y, pos.z};
    const Vec3 v1{pos.x, pos.y - ep, pos.z};
    const Vec3 v2{pos.x, pos.y, pos.z - ep};
    const float baseDistance = getDistance(pos);
    return Vec3{baseDistance - getDistance(v0), baseDistance - getDistance(v1),
                baseDistance - getDistance(v2)}
        .normalized();
  }

private:
  bool blocks[DEPTH][WIDTH]{};
};

Maze maze;

class Space {
public:
  Vec3 camPos{};
  float angle{};

  Vec3 renderPixel(int fragX, int fragY) const {
    Vec3 ray =
        Vec3{static_cast<float>(fragX * 2 - Canvas::WIDTH) / Canvas::HEIGHT,
             static_cast<float>(fragY * 2 - Canvas::HEIGHT) / Canvas::HEIGHT,
             -2.0f}
            .normalized()
            .rotateY(angle);
    const Vec3 light{-std::sin(angle), 0.0f, std::cos(angle)};
    Vec3 cur{camPos};
    for (int i = 0; i < 128; ++i) {
      const float d = getDistance(cur);
      if (d < 0.001f) {
        const auto normal = getNormal(cur);
        const float c = (0.5f + 0.5f * normal.dot(light)) *
                        (std::max(1.0f - (cur.sub(camPos)).length() / 4, 0.0f));
        return {c, c, c};
      }
      cur = cur.add(ray.mul(d));
    }
    return {0.0f, 0.0f, 0.0f};
  }

private:
  float getDistance(const Vec3 &pos) const {
    return maze.getDistance(pos);
#if 0
    auto p = pos.mod(6.0f).sub({3.0f, 3.0f, 3.0f});
    return p.length() - 1.0f;
#endif
  }

  Vec3 getNormal(const Vec3 &pos) const {
    const float ep = 0.0001f;
    const Vec3 v0{pos.x - ep, pos.y, pos.z};
    const Vec3 v1{pos.x, pos.y - ep, pos.z};
    const Vec3 v2{pos.x, pos.y, pos.z - ep};
    const float baseDistance = getDistance(pos);
    return Vec3{baseDistance - getDistance(v0), baseDistance - getDistance(v1),
                baseDistance - getDistance(v2)}
        .normalized();
  }
};

const float pi = 3.14159265258979f;

void update(Canvas &canvas) {
  static Vec3 pos{1.5f, 0.0f, 1.5f};
  static float angle = pi;
  static int x = 0;
  static int y = 0;
  static int vxc = 0;
  static int vyc = 0;
  const int sv = 1;
  if (InputManager::key(' '))
    vxc = vyc = 0;
  if (InputManager::key('a'))
    vxc = (vxc <= 0 ? -sv : 0);
  if (InputManager::key('d'))
    vxc = (vxc >= 0 ? sv : 0);
  if (InputManager::key('w'))
    vyc = (vyc >= 0 ? +sv : 0);
  if (InputManager::key('s'))
    vyc = (vyc <= 0 ? -sv : 0);
  float bfVel =
      (1.0f / 16) * ((vyc < 0 ? -1.0f : 0.0f) + (vyc > 0 ? 1.0f : 0.0f));
  float angleVel =
      (pi / 32) * ((vxc < 0 ? -1.0f : 0.0f) + (vxc > 0 ? 1.0f : 0.0f));
  angle += angleVel;
  pos.x += bfVel * cos(angle - 0.5f * pi);
  pos.z += bfVel * sin(angle - 0.5f * pi);
  float d = maze.getDistance(pos);
  const float r = 0.25f;
  if (d < r)
    pos = pos.add(maze.getNormal(pos).mul(r - d));
  canvas.clear();
  Space space;
  space.camPos.x = pos.x;
  space.camPos.z = pos.z;
  space.angle = angle;
  for (int j = 0; j < Canvas::HEIGHT; ++j) {
    for (int i = 0; i < Canvas::WIDTH; ++i) {
      const auto vc = space.renderPixel(i, j);
      canvas.pixel(i, j) = {static_cast<uint8_t>(vc.x * 255.9f),
                            static_cast<uint8_t>(vc.y * 255.9f),
                            static_cast<uint8_t>(vc.z * 255.9f)};
    }
  }
  canvas.pixel(static_cast<int>(pos.x + 10.0f),
               static_cast<int>(pos.z + 10.0f)) = {255, 0, 0};
}

int main() {
  std::srand(std::time(nullptr));
  maze.initialize();
#if 0
  printf("%d\n", std::rand() % 4);
  printf("%d\n", std::rand() % 4);
  printf("%d\n", std::rand() % 4);
  exit(1);
#endif
  Terminal::initialize();
  auto canvas = std::make_unique<Canvas>();
  int64_t timeMs = getCurrentTimeMs();
  while (!Terminal::exitRequested) {
    InputManager::update();
    update(*canvas);
    timeMs += 100;
    sleepTillMs(timeMs);
    canvas->printToTerminal();
  }
  Terminal::finalize();
}
