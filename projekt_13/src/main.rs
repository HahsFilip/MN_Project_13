
fn internal_detection(i: i32, j: i32, a: i32, b:i32, c:i32, d:i32) -> i32 {
    if (j < d-1 && i < b) || (j < d && i >= a - b){
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
                }
                if domain[i][j] != -2{
                    tmp = tmp + x[i][j]*(1.+4.*alpha+alpha*gamma);

                }else{
                    tmp = tmp + x[i][j]*(1.+4.*alpha);
                }
            }
            result[i][j] = tmp;
        }
        
    }
    result
} 
fn main() {

    const A: usize = 13;
    const B: usize = 5;
    const C: usize = 10;
    const D: usize = 8;
    let a_int: i32 = A.try_into().unwrap();
    let b_int: i32 = B.try_into().unwrap();
    let c_int: i32 = C.try_into().unwrap();
    let d_int: i32 = D.try_into().unwrap();

    let mut control_array = vec![1.0; 2*A-2*B];

    let mut domain_spec= vec![vec![0; C]; A];
    
    let mut u= vec![vec![1.0; C]; A];
    let mut a_matrix = vec![vec![0.0; C]; A];
    let mut test = vec![vec![0.0; C]; A];
    for i in 0..A{
        for j in 0..C{
            domain_spec[i][j] = internal_detection(i.try_into().unwrap(), j.try_into().unwrap(), a_int,b_int,c_int,d_int);
            
        }
        println!("{:?}",domain_spec[i]);
    }

    test = multiply_by_a_matrix(&mut u,&mut domain_spec,1.1,1.1);
    for i in 0..A{

        println!("{:?}",test[i]);
    }


}
