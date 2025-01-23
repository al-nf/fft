use num_complex::Complex;
use std::f64::consts::PI;

fn main() {
    let mut f: Vec<Complex<f64>> = Vec::new();
    for i in 1..5 {
        f.push(Complex::new(i as f64, 0.0))
    }
    for _ in 0..4 {
        f.push(Complex::new(0.0, 0.0));
    }

    let transformed = fast_fourier_transform(&f, 1.0);

    for num in transformed {
        println!("({:.5})", num);
    }
}

// Performs a DFT on a function <f> with the normalization factor <d>
fn fast_fourier_transform(f: &Vec<Complex<f64>>, d: f64) -> Vec<Complex<f64>> {
    // Calculates the base-2 logarithm of a decimal number <n> to find the number of bits needed to
    // represent it
    let log2 = |n: i32| -> i32 {
        let mut i = 0;
        let mut k = n;
        while k > 1 {
            k >>= 1;
            i += 1;
        }
        return i;
    };

    // Bit-reverses the number at index <n> for a DFT of size <k>
    let reverse = |n: i32, k: i32| {
        let mut result: i32 = 0;
        let bits = log2(k);

        for i in 0..bits {
            if n & (1 << i) != 0 {
                result |= 1 << (bits - 1 - i);
            }
        }
        return result;
    };

    // Cooley-Tukey
    let transform = |f: &Vec<Complex<f64>>| -> Vec<Complex<f64>> {
        let mut out = f.clone();
        let k = f.len();

        for i in 0..k {
            let j = reverse(i as i32, k as i32);
            if (i as i32) < j {
                out.swap(i, j as usize);
            }
        }
        
        let mut w: Vec<Complex<f64>> = Vec::with_capacity(k/2);
        
        for i in 0..k/2 {
            w.push(Complex::from_polar(1.0, -2.0 * PI * i as f64 / k as f64));
        }

        let mut len = 2;
        while len <= k {
            let half_len = len / 2;
            for i in (0..k).step_by(len) {
                for j in 0..half_len {
                    let u = out[i + j];
                    let v = w[k / len * j] * out[i + j + half_len];
                    out[i + j] = u + v;
                    out[i + j + half_len] = u - v;
                }
            }
            len *= 2;
        }
        out
    };

    // Applies the normalization factor <d>
    let mut transformed = transform(f);
    for i in 0..transformed.len() {
        transformed[i] *= d;
    }
    transformed
}
