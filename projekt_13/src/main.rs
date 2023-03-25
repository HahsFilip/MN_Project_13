
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
fn main() {

    const A: usize = 13;
    const B: usize = 5;
    const C: usize = 10;
    const D: usize = 8;
    let A_int: i32 = A.try_into().unwrap();
    let B_int: i32 = B.try_into().unwrap();
    let C_int: i32 = C.try_into().unwrap();
    let D_int: i32 = D.try_into().unwrap();

    let mut domain_spec: [[i32;C];A] = [[0;C];A];
    for i in 0..A{
        for j in 0..C{
            domain_spec[i][j] = internal_detection(i.try_into().unwrap(), j.try_into().unwrap(), A_int,B_int,C_int,D_int);
            
        }
        println!("{:?}",domain_spec[i]);
    }

}
