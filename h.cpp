#include <iostream>
#include <vector>
#include <cmath>

class Complex {
private:
    double re_;
    double im_;
    friend Complex operator+(const Complex&, const Complex&);
    friend Complex operator-(const Complex&, const Complex&);
    friend Complex operator*(const Complex&, const Complex&);
public:
    explicit Complex(double r = 0.0, double i = 0.0): re_(r), im_(i) {};
    Complex(const Complex& c) = default;
 
    Complex& operator=(const Complex& c) = default;
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
 
 
    Complex operator-() const {
        return Complex(-re_, -im_);
    }
 
    Complex operator+() const {
        return {Complex(*this)};
    }
 
 
    bool operator==(const Complex& c) const {
        return (im_ == c.im_) && (re_ == c.re_);
    }
 
    bool operator!=(const Complex& c) const {
        return !(*this == c);
    }
 
    double GetRe() const {return re_;};
    double GetIm() const {return im_;};
};
 
Complex operator+(const Complex& a, const Complex& b) {
    return Complex(a.re_ + b.re_, a.im_ + b.im_);
}
 
Complex operator-(const Complex& a, const Complex& b) {
    return Complex(a.re_ - b.re_, a.im_ - b.im_);
}
 
Complex operator*(const Complex& a, const Complex& b) {
    return Complex(a.re_ * b.re_ - a.im_ * b.im_, a.re_ * b.im_ + a.im_ * b.re_);
}
 
template <typename U>
U Abs(U number) {
  return number > 0 ? number : -number;
}
 
double Abs(Complex number) {
    return std::sqrt(number.GetRe() * number.GetRe() + number.GetIm() * number.GetIm());
}

const double kPI = acos(-1.0);

void FFT(std::vector<Complex>& poly, int invert = 0) {
    int n = static_cast<int>(poly.size());
    for(int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for(; j&bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(poly[i], poly[j]);
        }
    }
    for(int len = 2; len <= n; len <<= 1) {
        double ang = 2 * kPI / len * (invert ? 1 : -1);
        Complex wlen(cos(ang), sin(ang));
        int len_divided = len >> 1;
        for(int i = 0; i < n; i += len) {
            Complex w(1);
            for(int j = 0; j < len_divided; ++j) {
                Complex a = poly[i + j];
                Complex b = poly[i + j + len_divided] * w;
                poly[i + j] = a + b;
                poly[i + j + len_divided] = a - b;
                w *= wlen;
            }
        }
    }
    if (invert) {
        Complex div(1./static_cast<double>(n));
        for(Complex& i: poly) {
            i *= div;
        }
    }
}

std::vector<int> Multiply(std::vector<int>& a, std::vector<int>& b) {
    std::vector<Complex> polynom_a(a.begin(), a.end());
    std::vector<Complex> polynom_b(b.begin(), b.end());
    int len = 1;
    while(len < static_cast<int>(a.size() + b.size())) {
        len <<= 1;
    }

    polynom_a.resize(len);
    polynom_b.resize(len);
    FFT(polynom_a);
    FFT(polynom_b);

    for(int i = 0; i < len; ++i) {
        polynom_a[i] *= polynom_b[i];
    }

    FFT(polynom_a, 1);

    std::vector<int> final_polynom(len);
    for(int i = 0; i < len; ++i) {
        final_polynom[i] = static_cast<int>(round(polynom_a[i].GetRe()));
    }
    return final_polynom;
}

int main() {
    int n = 0;
    std::cin >> n;
    std::vector<int> a(n + 1);
    for(int i = 0; i <= n; ++i) {
        scanf("%d", &a[i]);
    }
    int m = 0;
    std::cin >> m;
    std::vector<int> b(m + 1);
    for(int i = 0; i <= m; ++i) {
        scanf("%d", &b[i]);
    }
    std::vector<int> product = Multiply(a, b);
    std::cout << n + m << " ";
    for(int i = 0; i < n + m + 1; ++i) {
        printf("%d ", product[i]);
    }
    return 0;
}