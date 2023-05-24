// Meisell-Lehmer
#include <cmath>
#include <cstdio>
#include <iostream>

class Constant {
 public:
  Constant() = default;
  
  static const size_t kN = 5e6 + 2;
  static const size_t kM = 7;
  static const size_t kPM = 2 * 3 * 5 * 7 * 11 * 13 * 17;
};

class PrimeNumberHelper {
 public:
   PrimeNumberHelper() = default;

   static bool np[Constant::kN];
   static size_t prime[Constant::kN];
   static size_t pi[Constant::kN];
   static size_t phi[Constant::kPM + 1][Constant::kM + 1];
   static size_t sz[Constant::kM + 1];
};

// const size_t Constant::kN;
// bool prime_number_helper.np[Constant::kN];
// size_t prime_number_helper.prime[Constant::kN];
// size_t prime_number_helper.pi[Constant::kN];
// const size_t Constant::kM;
// const size_t Constant::kPM;
// size_t prime_number_helper.phi[Constant::kPM + 1][Constant::kM + 1];
// size_t prime_number_helper.sz[Constant::kM + 1];

size_t GeneratePrimes(PrimeNumberHelper& prime_number_helper, Constant& constant) {
  size_t cnt = 0;
  prime_number_helper.np[0] = prime_number_helper.np[1] = true;
  prime_number_helper.pi[0] = prime_number_helper.pi[1] = 0;
  for (size_t i = 2; i < constant.kN; ++i) {
    if (!prime_number_helper.np[i]) {
      prime_number_helper.prime[++cnt] = i;
    }
    prime_number_helper.pi[i] = cnt;

    for (size_t j = 1; j <= cnt && i * prime_number_helper.prime[j] < constant.kN;
         ++j) {
      prime_number_helper.np[i * prime_number_helper.prime[j]] = true;
      if (i % prime_number_helper.prime[j] == 0) {
        break;
      }
    }
  }
  return cnt;
}

size_t SqrtUniversal(size_t x, size_t degree) {
  size_t r = 0;
  if (degree == 2) {
    r = static_cast<size_t>(sqrt(static_cast<double>(x) - 0.1));
    while (r * r <= x) {
      ++r;
    }
  } else {
    r = static_cast<size_t>(cbrt(static_cast<double>(x) - 0.1));
    while (r * r * r <= x) {
      ++r;
    }
  }
  return r - 1;
}

size_t CountNumbersWithKPrimeFactors(size_t x, size_t s, PrimeNumberHelper& prime_number_helper, Constant& constant) {
  if (s == 0) {
    return x;
  }
  if (s <= constant.kM) {
    return prime_number_helper.phi[x % prime_number_helper.sz[s]][s] +
           (x / prime_number_helper.sz[s]) * prime_number_helper.phi[prime_number_helper.sz[s]][s];
  }
  if (x <= prime_number_helper.prime[s] * prime_number_helper.prime[s]) {
    return prime_number_helper.pi[x] - s + 1;
  }
  if (x <= prime_number_helper.prime[s] * prime_number_helper.prime[s] * prime_number_helper.prime[s] &&
      x < constant.kN) {
    size_t s2x = prime_number_helper.pi[SqrtUniversal(x, 2)];
    size_t ans = prime_number_helper.pi[x] - (s2x + s - 2) * (s2x - s + 1) / 2;
    for (size_t i = s + 1; i <= s2x; ++i) {
      ans += prime_number_helper.pi[x / prime_number_helper.prime[i]];
    }
    return ans;
  }

  return CountNumbersWithKPrimeFactors(x, s - 1, prime_number_helper, constant) - CountNumbersWithKPrimeFactors(x / prime_number_helper.prime[s], s - 1, prime_number_helper, constant);
}

size_t CountPrimesPiFunc(size_t x, PrimeNumberHelper& prime_number_helper, Constant& constant) {
  if (x < constant.kN) {
    return prime_number_helper.pi[x];
  }
  size_t ans = CountNumbersWithKPrimeFactors(x, prime_number_helper.pi[SqrtUniversal(x, 3)], prime_number_helper, constant) +
                prime_number_helper.pi[SqrtUniversal(x, 3)] - 1;
  for (size_t i = prime_number_helper.pi[SqrtUniversal(x, 3)] + 1,
               ed = prime_number_helper.pi[SqrtUniversal(x, 2)];
       i <= ed; ++i) {
    ans -= CountPrimesPiFunc(x / prime_number_helper.prime[i], prime_number_helper, constant) - i + 1;
  }
  return ans;
}

size_t CountPrimes(size_t x, PrimeNumberHelper& prime_number_helper, Constant& constant) {
  if (x < constant.kN) {
    return prime_number_helper.pi[x];
  }
  size_t a = CountPrimes(SqrtUniversal(SqrtUniversal(x, 2), 2), prime_number_helper, constant);
  size_t b = CountPrimes(SqrtUniversal(x, 2), prime_number_helper, constant);
  size_t c = CountPrimes(SqrtUniversal(x, 3), prime_number_helper, constant);
  size_t sum =
      CountNumbersWithKPrimeFactors(x, a, prime_number_helper, constant) + static_cast<size_t>(b + a - 2) * (b - a + 1) / 2;
  for (size_t i = a + 1; i <= b; i++) {
    size_t w = x / prime_number_helper.prime[i];
    sum -= CountPrimes(w, prime_number_helper, constant);
    if (i > c) {
      continue;
    }
    size_t lim = CountPrimes(SqrtUniversal(w, 2), prime_number_helper, constant);
    for (size_t j = i; j <= lim; j++) {
      sum -= CountPrimes(w / prime_number_helper.prime[j], prime_number_helper, constant) - (j - 1);
    }
  }
  return sum;
}

int main() {

  size_t n = 0;
  std::cin >> n;

  Constant constant;
  PrimeNumberHelper prime_number_helper;

  GeneratePrimes(prime_number_helper, constant);

  prime_number_helper.sz[0] = 1;
  for (size_t i = 0; i <= constant.kPM; ++i) {
    prime_number_helper.phi[i][0] = i;
  }
  for (size_t i = 1; i <= constant.kM; ++i) {
    prime_number_helper.sz[i] = prime_number_helper.prime[i] * prime_number_helper.sz[i - 1];
    for (size_t j = 1; j <= constant.kPM; ++j) {
      prime_number_helper.phi[j][i] = prime_number_helper.phi[j][i - 1] -
                            prime_number_helper.phi[j / prime_number_helper.prime[i]][i - 1];
    }
  }

  std::cout << CountPrimes(n, prime_number_helper, constant);
  return 0;
}
