extern crate nalgebra as na;
use na::{SVector, SMatrix};


fn kalman_filter<const Z: usize, const X: usize>(z: &SVector<f64, Z>, r: &SMatrix<f64, Z, Z>,
                                                 x0: &SVector<f64, X>, p0: &SMatrix<f64, X,X>,
                                                h:&SMatrix<f64, Z,X>) -> (SVector<f64, X>, SMatrix<f64, X,X>)
{
   let y = h * x0;
   let residual = z - y;
   let p_inv = p0.try_inverse().unwrap();
   let s_inv = (h * p_inv * h.transpose() + r).try_inverse().unwrap();
   let k = p_inv * h.transpose() * s_inv;
   let x = x0 + k * residual;
   let p = (SMatrix::<f64, X,X>::identity() - k * h) * p_inv *
                                               (SMatrix::<f64, X,X>::identity() - k * h).transpose() + k * r * k.transpose();
    (x, p)
}

type x_t = SVector<f64, 6>;
type p_t = SMatrix<f64, 6,6>;

fn main() {

    let x0 = x_t::zeros();
    let p0 = p_t::from_diagonal_element(500.0);

    
    
    println!("Hello, world!");
}

