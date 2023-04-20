//Meisell-Lehmer
#include <cstdio>
#include <cmath>
#include <iostream>

const int64_t kN = 5e6 + 2;
bool np[kN];
int64_t prime[kN], pi[kN];
int64_t GetPrime()
{
    int64_t cnt = 0;
    np[0] = np[1] = true;
    pi[0] = pi[1] = 0;
    for(int64_t i = 2; i< kN; ++i) {
        if(!np[i])  {
            prime[++cnt] = i;
        }
        pi[i] = cnt;

        for(int64_t j = 1; j <= cnt && i * prime[j] < kN; ++j) {
            np[i * prime[j]] = true;
            if(i % prime[j] == 0) {
                break;
            }
        }
    }
    return cnt;
}

const int64_t kM= 7;
const int64_t kPM= 2 * 3 * 5 * 7 * 11 * 13 * 17;
int64_t phi[kPM+ 1][kM+ 1];
int64_t sz[kM+ 1];

void Init() {
    GetPrime();
    sz[0] = 1;
    for(int64_t i = 0; i <= kPM; ++i) {
        phi[i][0] = i;
    }
    for(int64_t i = 1; i <= kM; ++i) {
        sz[i] = prime[i] * sz[i - 1];
        for(int64_t j = 1; j <= kPM; ++j) {
            phi[j][i] = phi[j][i - 1] - phi[j / prime[i]][i - 1];
        }
    }
}

int64_t Sqrt2(int64_t x)
{
    auto r = static_cast<int64_t>(sqrt(static_cast<double>(x) - 0.1));
    while(r * r <= x) {
        ++r;
    }
    return r - 1;
}

int64_t Sqrt3(int64_t x) {
    auto r = static_cast<int64_t>(cbrt(static_cast<double>(x) - 0.1));
    while(r * r * r <= x) {
        ++r;
    }
    return r - 1;
}

int64_t GetPhi(int64_t x, int64_t s) {
    if(s == 0) {
        return x;
    } 
    if(s <= kM) {
        return phi[x % sz[s]][s] + (x / sz[s]) * phi[sz[s]][s];
    }
    if(x <= prime[s]*prime[s]) {
        return pi[x] - s + 1;
    }
    if(x <= prime[s]*prime[s]*prime[s] && x < kN) {
        int64_t s2x = pi[Sqrt2(x)];
        int64_t ans = pi[x] - (s2x +s - 2) * (s2x - s +1) / 2;
        for(int64_t i = s + 1; i <= s2x; ++i) {
            ans += pi[x / prime[i]];
        }
        return ans;
    }

    return GetPhi(x, s - 1) - GetPhi(x / prime[s], s - 1);
}

int64_t GetPi(int64_t x) {
    if(x < kN) {
        return pi[x];
    }
    int64_t ans = GetPhi(x, pi[Sqrt3(x)]) + pi[Sqrt3(x)] - 1;
    for(int64_t i = pi[Sqrt3(x)] + 1, ed = pi[Sqrt2(x)]; i <= ed; ++i) {
        ans -= GetPi(x / prime[i]) - i + 1;
    }
    return ans;
}

int64_t LehmerPi(int64_t x) {
    if(x < kN) {
        return pi[x];
    }
    int64_t a = LehmerPi(Sqrt2(Sqrt2(x)));
    int64_t b = LehmerPi(Sqrt2(x));
    int64_t c = LehmerPi(Sqrt3(x));
    int64_t sum = GetPhi(x, a) + static_cast<int64_t>(b +a - 2) * (b-a +1) / 2;
    for (int64_t i =a +1; i <= b; i++) {
        int64_t w = x / prime[i];
        sum -= LehmerPi(w);
        if (i > c) {
            continue;
        } 
        int64_t lim = LehmerPi(Sqrt2(w));
        for (int64_t j = i; j <= lim; j++) {
            sum -= LehmerPi(w / prime[j]) - (j - 1);
        }
    }
    return sum;
}

int main()
{
    Init();
    int64_t n = 0;
    std::cin >> n;
    std::cout << LehmerPi(n);
    return 0;
}