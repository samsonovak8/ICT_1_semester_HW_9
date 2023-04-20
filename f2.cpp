// Meisell-Lehmer
#include <cmath>
#include <cstdio>
#include <iostream>

class Constant {
 public:
  static const int64_t kN = 5e6 + 2;
  static bool np[kN];
  static int64_t prime[kN];
  static int64_t pi[kN];
  static const int64_t kM = 7;
  static const int64_t kPM = 2 * 3 * 5 * 7 * 11 * 13 * 17;
  static int64_t phi[kPM + 1][kM + 1];
  static int64_t sz[kM + 1];
};

const int64_t Constant::kN;
bool Constant::np[Constant::kN];
int64_t Constant::prime[Constant::kN];
int64_t Constant::pi[Constant::kN];
const int64_t Constant::kM;
const int64_t Constant::kPM;
int64_t Constant::phi[Constant::kPM + 1][Constant::kM + 1];
int64_t Constant::sz[Constant::kM + 1];

int64_t GetPrime() {
  int64_t cnt = 0;
  Constant::np[0] = Constant::np[1] = true;
  Constant::pi[0] = Constant::pi[1] = 0;
  for (int64_t i = 2; i < Constant::kN; ++i) {
    if (!Constant::np[i]) {
      Constant::prime[++cnt] = i;
    }
    Constant::pi[i] = cnt;

    for (int64_t j = 1; j <= cnt && i * Constant::prime[j] < Constant::kN;
         ++j) {
      Constant::np[i * Constant::prime[j]] = true;
      if (i % Constant::prime[j] == 0) {
        break;
      }
    }
  }
  return cnt;
}

void Init() {
  GetPrime();
  Constant::sz[0] = 1;
  for (int64_t i = 0; i <= Constant::kPM; ++i) {
    Constant::phi[i][0] = i;
  }
  for (int64_t i = 1; i <= Constant::kM; ++i) {
    Constant::sz[i] = Constant::prime[i] * Constant::sz[i - 1];
    for (int64_t j = 1; j <= Constant::kPM; ++j) {
      Constant::phi[j][i] = Constant::phi[j][i - 1] -
                            Constant::phi[j / Constant::prime[i]][i - 1];
    }
  }
}

int64_t SqrtUniversal(int64_t x, int64_t degree) {
  int64_t r = 0;
  if (degree == 2) {
    r = static_cast<int64_t>(sqrt(static_cast<double>(x) - 0.1));
    while (r * r <= x) {
      ++r;
    }
  } else {
    r = static_cast<int64_t>(cbrt(static_cast<double>(x) - 0.1));
    while (r * r * r <= x) {
      ++r;
    }
  }
  return r - 1;
}

int64_t GetPhi(int64_t x, int64_t s) {
  if (s == 0) {
    return x;
  }
  if (s <= Constant::kM) {
    return Constant::phi[x % Constant::sz[s]][s] +
           (x / Constant::sz[s]) * Constant::phi[Constant::sz[s]][s];
  }
  if (x <= Constant::prime[s] * Constant::prime[s]) {
    return Constant::pi[x] - s + 1;
  }
  if (x <= Constant::prime[s] * Constant::prime[s] * Constant::prime[s] &&
      x < Constant::kN) {
    int64_t s2x = Constant::pi[SqrtUniversal(x, 2)];
    int64_t ans = Constant::pi[x] - (s2x + s - 2) * (s2x - s + 1) / 2;
    for (int64_t i = s + 1; i <= s2x; ++i) {
      ans += Constant::pi[x / Constant::prime[i]];
    }
    return ans;
  }

  return GetPhi(x, s - 1) - GetPhi(x / Constant::prime[s], s - 1);
}

int64_t GetPi(int64_t x) {
  if (x < Constant::kN) {
    return Constant::pi[x];
  }
  int64_t ans = GetPhi(x, Constant::pi[SqrtUniversal(x, 3)]) +
                Constant::pi[SqrtUniversal(x, 3)] - 1;
  for (int64_t i = Constant::pi[SqrtUniversal(x, 3)] + 1,
               ed = Constant::pi[SqrtUniversal(x, 2)];
       i <= ed; ++i) {
    ans -= GetPi(x / Constant::prime[i]) - i + 1;
  }
  return ans;
}

int64_t LehmerPi(int64_t x) {
  Init();
  if (x < Constant::kN) {
    return Constant::pi[x];
  }
  int64_t a = LehmerPi(SqrtUniversal(SqrtUniversal(x, 2), 2));
  int64_t b = LehmerPi(SqrtUniversal(x, 2));
  int64_t c = LehmerPi(SqrtUniversal(x, 3));
  int64_t sum =
      GetPhi(x, a) + static_cast<int64_t>(b + a - 2) * (b - a + 1) / 2;
  for (int64_t i = a + 1; i <= b; i++) {
    int64_t w = x / Constant::prime[i];
    sum -= LehmerPi(w);
    if (i > c) {
      continue;
    }
    int64_t lim = LehmerPi(SqrtUniversal(w, 2));
    for (int64_t j = i; j <= lim; j++) {
      sum -= LehmerPi(w / Constant::prime[j]) - (j - 1);
    }
  }
  return sum;
}

int main() {
  int64_t n = 0;
  std::cin >> n;
  std::cout << LehmerPi(n);
  return 0;
}
