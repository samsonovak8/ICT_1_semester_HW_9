#include <cmath>
#include <iostream>
#include <vector>

const double kPI = acos(-1.0);

class Complex {
 public:
  explicit Complex(double r = 0.0, double i = 0.0) : re_(r), im_(i){};
  ~Complex() noexcept = default;

  Complex(const Complex& c) = default;

  Complex& operator=(const Complex& c) = default;

  Complex(Complex&& c) noexcept : re_(c.re_), im_(c.im_) {
    c.re_ = 0;
    c.im_ = 0;
  }

  Complex& operator=(Complex&& c) = default;

  Complex& operator+=(const Complex& c) {
    re_ += c.re_;
    im_ += c.im_;
    return *this;
  }

  Complex& operator-=(const Complex& c) {
    re_ -= c.re_;
    im_ -= c.im_;
    return *this;
  }

  Complex& operator*=(const Complex& c) {
    double re_t = re_;
    double im_t = im_;
    re_ = re_t * c.re_ - im_t * c.im_;
    im_ = re_t * c.im_ + im_t * c.re_;
    return *this;
  }

  Complex operator-() const { return {-re_, -im_}; }

  Complex operator+() const { return {Complex(*this)}; }

  bool operator==(const Complex& c) const {
    return (im_ == c.im_) && (re_ == c.re_);
  }

  bool operator!=(const Complex& c) const { return !(*this == c); }

  double GetRe() const { return re_; };
  double GetIm() const { return im_; };

  Complex operator+(const Complex& b) {
    Complex tmp(*this);
    return tmp += b;
  }

  Complex operator-(const Complex& b) {
    Complex tmp(*this);
    return tmp -= b;
  }

  Complex operator*(const Complex& b) {
    Complex tmp(*this);
    return tmp *= b;
  }

  friend std::istream& operator>>(std::istream& is, Complex& c) {
    is >> c.re_ >> c.im_;
    return is;
  }

  friend std::ostream& operator<<(std::ostream& os, const Complex& c) {
    os << c.re_ << " " << c.im_ << " ";
    return os;
  }

 private:
  double re_;
  double im_;
};

template <typename T>
T Abs(T number) {
  return number > 0 ? number : -number;
}

double Abs(Complex number) {
  return std::sqrt(number.GetRe() * number.GetRe() +
                   number.GetIm() * number.GetIm());
}

void FFT(std::vector<Complex>& poly, int invert = 0) {
  int n = static_cast<int>(poly.size());
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      std::swap(poly[i], poly[j]);
    }
  }
  for (int len = 2; len <= n; len <<= 1) {
    double ang = 2 * kPI / len * (invert ? 1 : -1);
    Complex wlen(cos(ang), sin(ang));
    int len_divided = len >> 1;
    for (int i = 0; i < n; i += len) {
      Complex w(1);
      for (int j = 0; j < len_divided; ++j) {
        Complex a = poly[i + j];
        Complex b = poly[i + j + len_divided] * w;
        poly[i + j] = a + b;
        poly[i + j + len_divided] = a - b;
        w *= wlen;
      }
    }
  }
  if (invert) {
    Complex div(1. / static_cast<double>(n));
    for (Complex& i : poly) {
      i *= div;
    }
  }
}

class Polynom {
 public:
  Polynom() {}

  Polynom(const std::vector<Complex>& coeffs, size_t degree)
      : coefficients(coeffs), degree_(degree) {}

  Polynom operator*(const Polynom& other) {
    // std::vector<double> resultCoeffs(coefficients.size() +
    // other.coefficients.size() - 1, 0.0);

    // for (size_t i = 0; i < coefficients.size(); ++i) {
    //     for (size_t j = 0; j < other.coefficients.size(); ++j) {
    //         resultCoeffs[i + j] += coefficients[i] * other.coefficients[j];
    //     }
    // }

    // return Polynom(resultCoeffs, resultCoeffs.size());
    std::vector<Complex> polynom_a(coefficients.begin(), coefficients.end());
    std::vector<Complex> polynom_b(other.coefficients.begin(),
                                   other.coefficients.end());
    int len = 1;
    int up_border = static_cast<int>(degree_ + other.degree_);
    while (len < up_border) {
      len <<= 1;
    }

    polynom_a.resize(len);
    polynom_b.resize(len);
    FFT(polynom_a);
    FFT(polynom_b);

    for (int i = 0; i < len; ++i) {
      polynom_a[i] *= polynom_b[i];
    }

    FFT(polynom_a, 1);

    std::vector<int> final_polynom(len);
    for (int i = 0; i < len; ++i) {
      final_polynom[i] = static_cast<int>(round(polynom_a[i].GetRe()));
    }
    return {final_polynom, final_polynom.size()};
  }

  friend std::istream& operator>>(std::istream& is, Polynom& p) {
    size_t degree;
    is >> degree;
    p.degree_ = degree;
    p.coefficients.resize(degree + 1);
    for (size_t i = 0; i <= degree; ++i) {
      is >> p.coefficients[i];
    }

    return is;
  }

  size_t GetDegree() { return degree_; }

  friend std::ostream& operator<<(std::ostream& os, const Polynom& p) {
    Complex zero(0.0, 0.0);
    for (int i = p.coefficients.size() - 1; i >= 0; --i) {
      if (p.coefficients[i] != zero) {
        os << p.coefficients[i] << " ";
      }
    }
    return os;
  }

 private:
  friend class Complex;
  std::vector<Complex> coefficients;
  size_t degree_;
};

int main() {
  Polynom a;
  std::cin >> a;
  Polynom b;
  std::cin >> b;
  size_t n = a.GetDegree();
  size_t m = b.GetDegree();
  Polynom product = a * b;
  std::cout << n + m << " ";
  std::cout << product;
  return 0;
}
