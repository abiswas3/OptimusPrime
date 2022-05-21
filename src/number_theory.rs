// Just for scalars right now: Packed needs some optimistions 



/// Euclids algorithm (recursive from https://www.wikiwand.com/en/Euclidean_algorithm)
pub fn gcd(a: i64, b: i64)->(i64, i64, i64){

    if b == 0 {
        return (a, 1, 0)
    }

    let n = a / b;
    let c =  a % b;
    let r = gcd(b, c);
    return (r.0, r.2, r.1 - r.2 * n)

}

#[test]
fn test_gcd() {
    assert_eq!(gcd(12, 16), (4, -1, 1));
}

/// Inverse of `k` in the *Zp* field defined by `prime`.
pub fn mod_inverse(k: i64, prime: i64) -> i64 {
    let k2 = k % prime;
    let r = if k2 < 0 {
        -gcd(prime, -k2).2
    } else {
        gcd(prime, k2).2
    };
    (prime + r) % prime
}

#[test]
fn test_mod_inverse() {
    assert_eq!(mod_inverse(3, 7), 5);
}


/// `x` to the power of `e` in the *Zp* field defined by `prime`.
pub fn mod_pow(mut x: i64, mut e: u32, prime: i64) -> i64 {
    let mut acc = 1;
    while e > 0 {
        if e % 2 == 0 {
            // even
            // no-op
        } else {
            // odd
            acc = (acc * x) % prime;
        }
        x = (x * x) % prime; // waste one of these by having it here but code is simpler (tiny bit)
        e = e >> 1;
    }
    acc
}

#[test]
fn test_mod_pow() {
    assert_eq!(mod_pow(2, 0, 17), 1);
    assert_eq!(mod_pow(2, 3, 17), 8);
    assert_eq!(mod_pow(2, 6, 17), 13);

    assert_eq!(mod_pow(-3, 0, 17), 1);
    assert_eq!(mod_pow(-3, 1, 17), -3);
    assert_eq!(mod_pow(-3, 15, 17), -6);
}

pub fn mod_evaluate_polynomial(coefficients: &[i64], point: i64, prime: i64) -> i64 {
    // evaluate using Horner's rule
    //  - to combine with fold we consider the coefficients in reverse order
    // first element is the secret, last element is the coeffiecient with the highest deg
    let reversed_coefficients = coefficients.iter().rev();
    reversed_coefficients.fold(0, |partial, coef| (partial * point + coef).rem_euclid( prime))
}

#[test]
fn test_mod_evaluate_polynomial(){

    // p(x) = (2 + x + 2*x**2 + 3*x**3 ) % 11
    // let secret: i64 = 3;
    // THis is just for me to understand what's going on
    let prime = 11;
    let a = [2, 1, 2, 3];
    let point: i64 = 3;
    let sum = a.iter().rev().fold(0, |acc, coef| (acc * point + coef) % prime);
    assert_eq!(sum, 5);

    assert_eq!(mod_evaluate_polynomial(&a, point, prime), 5);

    // p(x) = (1 + 2*x + 3*x**2 + 4*x**3 + 5*x**4 + 6*x**5) % 17
    let poly = vec![1i64, 2, 3, 4, 5, 6];
    let point = 0;
    let prime = 17;
    assert_eq!(mod_evaluate_polynomial(&poly, point, prime), 1i64);

}

pub fn lagrange_interpolation_at_zero(points: &[i64], values: &[i64], prime: i64) -> i64 {
    assert_eq!(points.len(), values.len());
    // Lagrange interpolation for point 0
    let mut acc = 0i64;
    for i in 0..values.len() {
        let xi = points[i];
        let yi = values[i];
        let mut num = 1i64;
        let mut denum = 1i64;
        for j in 0..values.len() {
            if j != i {
                let xj = points[j];
                num = (num * xj) % prime;
                denum = (denum * (xj - xi)) % prime;
            }
        }
        acc = (acc + yi * num * mod_inverse(denum, prime)) % prime;
    }
    acc
}