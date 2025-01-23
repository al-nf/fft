#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
using namespace std;

// Calculates the base-2 logarithm of a decimal number <n> to find the number of bits needed to represent it
int log2(int n) {
    int i = 0;
    while (n > 1) {
        n >>= 1;
        i++;
    }
    return i;
}

// Bit-reverses the number at index <n> for a DFT of size <k>
int reverse(int n, int k) {
    int result = 0;
    int bits = log2(k);
    for (int i = 0; i < bits; i++) {
        if (n & (1 << i)) {
            result |= (1 << (bits - 1 - i));
        }
    }
    return result;
}

// Cooley-Tukey
void transform(vector<complex<double> >& f) {
    int k = f.size();

    for (int i = 0; i < k; i++) {
        int j = reverse(i, k);
        if (i < j) {
            swap(f[i], f[j]);
        }
    }

    vector<complex<double> > W(k / 2);
    for (int i = 0; i < k / 2; i++) {
        W[i] = polar(1.0, -2.0 * M_PI * i / k);
    }

    // FFT
    for (int len = 2; len <= k; len *= 2) {
        int halfLen = len / 2;
        for (int i = 0; i < k; i += len) {
            for (int j = 0; j < halfLen; j++) {
                complex<double> u = f[i + j];
                complex<double> v = W[k / len * j] * f[i + j + halfLen];
                f[i + j] = u + v;
                f[i + j + halfLen] = u - v;
            }
        }
    }
}

void FFT(vector<complex<double> >& f, double d) {
    transform(f);
    for (auto& val : f) {
        val *= d;
    }
}

int main() {
    vector<complex<double> > f;
    f.push_back(complex<double>(1, 0));
    f.push_back(complex<double>(2, 0));
    f.push_back(complex<double>(3, 0));
    f.push_back(complex<double>(4, 0));
    f.push_back(complex<double>(0, 0));
    f.push_back(complex<double>(0, 0));
    f.push_back(complex<double>(0, 0));
    f.push_back(complex<double>(0, 0));

    FFT(f, 1.0);

    for (const auto& val : f) {
        cout << val << '\n';
    }
}
