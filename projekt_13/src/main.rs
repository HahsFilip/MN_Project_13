use std::convert::TryFrom;

fn internal_detection(i: i32, j: i32, a: i32, b:i32, c:i32, d:i32) -> i32 {
    if (j < d-1 && i < b) || (j < d-1 && i >= a - b){
        -1 // outside
    }else if j == 0 && i >= b && i < a - b{
        i - b // get index of c on bottom of t
    }else if j == c-1{
       // println!("{}", i + (a-2*b));
        i + (a-2*b) // get index of c value
    }else{
        -2 // inside
    }

}
fn multiply_by_a_matrix(x: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>, mut alpha: f32, gamma: f32 ) -> Vec<Vec<f32>> {
    let mut tmp = 0.0;
    let mut result = vec![vec![0.0; x[0].len()]; x.len()];
    for i in 0..x.len(){
        for j in 0..x[0].len(){
            tmp = 0.0;
            if domain[i][j] != -1{
                if i != 0 && i < x.len() - 1 && j != 0 && j < x[0].len() - 1{
                    // internal area
                    if domain[i][j] == -2{
                        if domain[i+1][j] == -1{
                            tmp = tmp - alpha * x[i - 1][j];
                        }else{
                            tmp = tmp - alpha * x[i+1][j];
                        }
                        if domain[i-1][j] == -1{
                            tmp = tmp - alpha * x[i + 1][j];
                        }else{
                            tmp = tmp - alpha * x[i - 1][j];
                        }
                        if domain[i][j+1] == -1{
                            tmp = tmp - alpha * x[i][j-1];
                        }else{
                            tmp = tmp - alpha * x[i][j+1];
                        }
                        if domain[i][j-1] == -1{
                            tmp = tmp - alpha * x[i][j+1];
                        }else{
                            tmp = tmp - alpha * x[i][j-1];
                        }
                    
                    }
                }else{
                    if j == 0 {
                        if i > 0 && i < x.len() - 1 {
                            if domain[i+1][j] == -1{
                                tmp = tmp - alpha * x[i - 1][j];
                            }else{
                                tmp = tmp - alpha * x[i+1][j];
                            }
                            if domain[i-1][j] == -1{
                                tmp = tmp - alpha * x[i + 1][j];
                            }else{
                                tmp = tmp - alpha * x[i - 1][j];
                            }
                            tmp = tmp - 2.*alpha * x[i][j + 1];

                    }else{
                        panic!("here");

                    }


                        tmp = tmp + x[i][j]*(1.+4.*alpha+alpha*gamma);
                    }else if j == x[0].len()-1{
                        if (i > 0 && i < x.len()-1){
                            
                            if domain[i+1][j] == -1{
                                tmp = tmp - alpha * x[i - 1][j];
                            }else{
                                tmp = tmp - alpha * x[i + 1][j];
                            }
                            if domain[i-1][j] == -1{
                                tmp = tmp - alpha * x[i + 1][j];
                            }else{
                                tmp = tmp - alpha * x[i - 1][j];
                            }
                           // tmp = tmp - 2.*alpha * x[i][j - 1];
                        }else{
                            if i == 0{
                                tmp = tmp - 2.0*alpha*x[i+1][j] - alpha*x[i][j-1] ;
                                 
                        }else{
                            tmp = tmp - 2.0*alpha*x[i-1][j] - alpha*x[i][j-1];
                            //println!("{} {} {}",i,j, domain[i][j]); 
                        }
                        }
                        tmp = tmp + x[i][j]*(1.+4.*alpha+alpha*gamma);
                      
                    }else{
                        if i == 0{
                            tmp = tmp - 2.0*alpha*x[i+1][j] - alpha*x[i][j+1]-alpha*x[i][j-1] ;
                        }else{
                            tmp = tmp - 2.0*alpha*x[i-1][j] - alpha*x[i][j+1]-alpha*x[i][j-1];
                     
                        }
                    }
                    

                }
   
            }
            result[i][j] = tmp;
        }
        
    }
    return result;
}
fn compute_b (u: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>, c: &mut Vec<f32>, gamma: f32 ) -> Vec<Vec<f32>>{
    let mut result = vec![vec![0.0; u[0].len()]; u.len()];
    for i in 0..u.len(){
        for j in 0..u[0].len(){
            if domain[i][j] != -1{
                if domain[i][j] == -2{
                    result[i][j] = u[i][j];

                }else{
                    {
                        let index: usize = domain[i][j] as usize;
                    result[i][j] = u[i][j]+c[index]*gamma;
                    }
                }
            }
        }
    }
    return result;
}
fn subtrac_vec (x: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>, y: &mut Vec<Vec<f32>>)-> Vec<Vec<f32>>{// x - y
    let mut result = vec![vec![0.0; x[0].len()]; x.len()];
    for i in 0..x.len(){
        for j in 0..x[0].len(){
            if domain[i][j] != -1{
                result[i][j] = x[i][j] - y[i][j];
            }
        }
    }
    return result;
}
fn multiply_by_scalar_vec (x: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>, y: f32)-> Vec<Vec<f32>>{// x*y
    let mut result = vec![vec![0.0; x[0].len()]; x.len()];
    for i in 0..x.len(){
        for j in 0..x[0].len(){
            if domain[i][j] != -1{
                result[i][j] = x[i][j]*y;
            }
        }
    }
    return result;
}
fn scalar_product (x: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>,y: &mut Vec<Vec<f32>>)-> f32{
    let mut result: f32 = 0.0;
    for i in 0..x.len(){
        for j in 0..x[0].len(){
            if domain[i][j] != -1{
                result = result + x[i][j]*y[i][j];
            }
        }
    }
    return result;
}
fn scalar_product_itself (x: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>)-> f32{
    let mut result: f32 = 0.0;
    for i in 0..x.len(){
        for j in 0..x[0].len(){
            if domain[i][j] != -1{
                result = result + x[i][j]*x[i][j];
            }
        }
    }
    return result;
}
fn pretty_print_vec( x: &mut Vec<Vec<f32>>) {
    for i in 0..x.len(){
        println!("{:?}",x[i]);
    } 
}

fn main() {

    const A: usize = 20;
    const B: usize = 5;
    const C: usize = 15;
    const D: usize = 8;
    let a_int: i32 = A.try_into().unwrap();
    let b_int: i32 = B.try_into().unwrap();
    let c_int: i32 = C.try_into().unwrap();
    let d_int: i32 = D.try_into().unwrap();

    let mut control_array = vec![200.0; 2*A-2*B];

    let mut domain_spec= vec![vec![0; C]; A];
    
    let mut u= vec![vec![1.0; C]; A];
    let mut a_matrix = vec![vec![0.0; C]; A];
    let mut ax = vec![vec![0.0; C]; A];
    let mut r = vec![vec![0.0; C]; A];
    let mut b = vec![vec![0.0; C]; A];
    let mut z = vec![vec![0.0; C]; A];
    let mut p = vec![vec![0.0; C]; A];
    let mut tmp = vec![vec![0.0; C]; A];
    let mut gamma_sim_par = 0.01;
    let mut alpha_sim_par = 0.01;
    let mut delta_solve = 0.0;
    let mut alpha_solve = 0.0;
    let mut beta_solve = 0.0;
    let mut gamma_solve = 0.0;

    for i in 0..A{
        for j in 0..C{
            domain_spec[i][j] = internal_detection(i.try_into().unwrap(), j.try_into().unwrap(), a_int,b_int,c_int,d_int);
            
        }
       // println!("{:?}",domain_spec[i]);
    }
    
    ax = multiply_by_a_matrix(&mut u,&mut domain_spec,alpha_sim_par,gamma_sim_par);
   // pretty_print_vec(&mut ax);
   // println!("-------------------\n");
    b = compute_b(&mut u, &mut domain_spec, &mut control_array,gamma_sim_par);
    r = subtrac_vec(&mut b, &mut domain_spec, &mut ax);
    p = r.clone();
   // pretty_print_vec(&mut p);
    println!("-------------------\n");
    delta_solve = scalar_product_itself(&mut r, &mut domain_spec );
    let mut gamma_zero = delta_solve;
    //pretty_print_vec(&mut r);
    for n in 0..10{
      //  pretty_print_vec(&mut u);
        z = multiply_by_a_matrix(&mut p,&mut domain_spec,alpha_sim_par,gamma_sim_par);
        //pretty_print_vec(&mut z);
        alpha_solve = -delta_solve/scalar_product(&mut p, &mut domain_spec, &mut z);
        //println!("{}", delta_solve);
        tmp = multiply_by_scalar_vec(&mut u, &mut domain_spec, alpha_solve);
        //pretty_print_vec(&mut tmp);
        u = subtrac_vec(&mut u, &mut domain_spec, &mut tmp);
        
        tmp = multiply_by_scalar_vec(&mut z, &mut domain_spec, -alpha_solve);
        r = subtrac_vec(&mut r, &mut domain_spec, &mut tmp);
        gamma_solve = scalar_product_itself(&mut r, &mut domain_spec);
        println!("{}", gamma_solve/gamma_zero);

        beta_solve = gamma_solve/delta_solve;
        tmp = multiply_by_scalar_vec(&mut p, &mut domain_spec, -beta_solve);
        p = subtrac_vec(&mut r, &mut domain_spec, &mut tmp);
        delta_solve = gamma_solve;
        //pretty_print_vec(&mut u);
    }



}
