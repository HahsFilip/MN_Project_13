use rand::Rng;
extern crate sdl2;

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::rect::Rect;


use std::{thread, time};


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
fn multiply_by_a_matrix( x: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>, alpha: f32, gamma: f32 ) -> Vec<Vec<f32>> {
    let mut tmp:f32;
    let mut result = vec![vec![0.0; x[0].len()]; x.len()];
    for i in 1..x.len()-1{
        for j in 1..x[0].len()-1{
            tmp = 0.0;
            if domain[i][j] != -1{
                if domain[i][j] == -2{
                    tmp = tmp + (1.0+4.0*alpha)*x[i][j];
                }else{
                    tmp = tmp + (1.0+4.0*alpha+alpha*gamma)*x[i][j];
                }

                if domain[i][j+1] != -1{
                    tmp = tmp - alpha*x[i][j+1];
                }else{
                    tmp = tmp - alpha*x[i][j-1];
                }
                if domain[i][j-1] != -1{
                    tmp = tmp - alpha*x[i][j-1];
                }else{
                    tmp = tmp - alpha*x[i][j+1];
                }
                if domain[i-1][j] != -1{
                    tmp = tmp - alpha*x[i-1][j];
                }else{
                    tmp = tmp - alpha*x[i+1][j];
                }
                if domain[i+1][j] != -1{
                    tmp = tmp - alpha*x[i+1][j];
                }else{
                    tmp = tmp - alpha*x[i-1][j];
                }
            }
            result[i][j] = tmp;
        }

    }
    return result;
}
fn compute_b (u: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>, c:&mut Vec<f32>, gamma: f32, alpha: f32 ) -> Vec<Vec<f32>>{
    let mut result = vec![vec![0.0; u[0].len()]; u.len()];
    for i in 0..u.len(){
        for j in 0..u[0].len(){
            if domain[i][j] != -1{
                if domain[i][j] == -2{
                    result[i][j] = u[i][j];

                }else{
                    {
                        let index: usize = domain[i][j] as usize;
                    result[i][j] = u[i][j]+c[index]*gamma*alpha;
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
fn multiply_by_scalar_vec (x: &mut Vec<Vec<f32>>,domain: &mut
Vec<Vec<i32>>, y: f32)-> Vec<Vec<f32>>{// x*y
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
fn scalar_product (x: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>,y:
&mut Vec<Vec<f32>>)-> f32{
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
fn scalar_product_itself (x: &mut Vec<Vec<f32>>,domain: &mut
Vec<Vec<i32>>)-> f32{
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
        println!("{:.1?}",x[i]);
    }
}

fn pretty_print_int( x: &mut Vec<Vec<i32>>) {
    for i in 0..x.len(){
        println!("{:?}",x[i]);
    }
}
fn force_copy(x: &mut Vec<Vec<f32>> ) -> Vec<Vec<f32>> {
    let mut result = vec![vec![0.0; x[0].len()]; x.len()];
        for i in 0..x.len(){
        for j in 0..x[0].len(){

                result[i][j] = x[i][j];

        }
    }
    return result;
}

fn conjugate_gradiant(u_0:  Vec<Vec<f32>>,domain_spec: &mut Vec<Vec<i32>>,a_func: &dyn Fn( &mut Vec<Vec<f32>>, &mut Vec<Vec<i32>>, f32, f32) -> Vec<Vec<f32>>,b_func: &dyn Fn( &mut Vec<Vec<f32>>, &mut Vec<Vec<i32>>, f32, f32,f32) -> Vec<Vec<f32>>, f32: gamma_sim_par, f32: alpha_sim_par  ) -> Vec<Vec<f32>>{
    let mut ax :Vec<Vec<f32>>;
    let mut u :Vec<Vec<f32>>;
    let mut r :Vec<Vec<f32>>;
    let mut b :Vec<Vec<f32>>;
    let mut z :Vec<Vec<f32>>;
    let mut p :Vec<Vec<f32>>;
    let mut tmp :Vec<Vec<f32>>;
    let mut delta_solve:f32;
    let mut alpha_solve:f32;
    let mut beta_solve :f32;
    let mut gamma_solve:f32;
    u = u_0.clone();
    ax = a_func( u,domain_spec,alpha_sim_par,gamma_sim_par);
    // pretty_print_vec(&mut ax);
    // println!("-------------------\n");
     b = b_func(&mut u_0,  domain_spec, control_array,gamma_sim_par, alpha_sim_par);
     r = subtrac_vec( &mut b, domain_spec, &mut ax);
     p = r.clone();
     delta_solve = scalar_product_itself(&mut r, &mut domain_spec );
     let gamma_zero = delta_solve;
     for _n in 0..10{
        z = multiply_by_a_matrix(&mut p,&mut domain_spec,alpha_sim_par, gamma_sim_par);
        //pretty_print_vec(&mut z);
        alpha_solve = -delta_solve/scalar_product(&mut p, &mut domain_spec, &mut z);
        //println!("{}", delta_solve);
        tmp = multiply_by_scalar_vec(&mut p, &mut domain_spec, alpha_solve);
        //pretty_print_vec(&mut tmp);
        *u_0 = subtrac_vec(u_0, &mut domain_spec, &mut tmp);

        tmp = multiply_by_scalar_vec(&mut z, &mut domain_spec,-alpha_solve);
        r = subtrac_vec(&mut r, &mut domain_spec, &mut tmp);
        gamma_solve = scalar_product_itself(&mut r, &mut domain_spec);
        println!("{}", gamma_solve/gamma_zero);
        let distance = scalar_product_itself(&mut r, &mut domain_spec);

        beta_solve = gamma_solve/delta_solve;
        tmp = multiply_by_scalar_vec(&mut p, &mut domain_spec, beta_solve);
        p = subtrac_vec(&mut r, &mut domain_spec, &mut tmp);
        delta_solve = gamma_solve;
        //pretty_print_vec(&mut u);
    }
     u_0
}

fn main()-> Result<(), String> {
    let diffusivity = 100.0;
    let h = 1.0;
    let beta = 100.0;
    let dt = 0.10;

    let sdl_context = sdl2::init()?;
    let video_subsystem = sdl_context.video()?;

    let window = video_subsystem
        .window("rust-sdl2 demo: Video", 800, 600)
        .position_centered()
        .opengl()
        .build()
        .map_err(|e| e.to_string())?;

    let mut canvas = window.into_canvas().build().map_err(|e| e.to_string())?;
    let mut event_pump = sdl_context.event_pump()?;

    canvas.set_draw_color(Color::RGB(0, 0, 0));


    canvas.present();


    const A: usize = 200;
    const B: usize = 50;
    const C: usize = 150;
    const D: usize = 80;
    let a_int: i32 = A.try_into().unwrap();
    let b_int: i32 = B.try_into().unwrap();
    let c_int: i32 = C.try_into().unwrap();
    let d_int: i32 = D.try_into().unwrap();

    let mut control_array = vec![100.0; 2*A-2*B];

    let mut domain_spec= vec![vec![-1; C+2]; A+2];
    let mut rng = rand::thread_rng();

    let mut u= vec![vec![1.1; C+2]; A+2];

    let gamma_sim_par = 2.0*h*beta/diffusivity;
    let alpha_sim_par = diffusivity*dt/(h*h);


    for i in 0..A{
        for j in 0..C{
            domain_spec[i+1][j+1] =
internal_detection(i.try_into().unwrap(), j.try_into().unwrap(),
a_int,b_int,c_int,d_int);
            u[i+1][j+1] = rng.gen::<f32>()*50.1;
        }
       // println!("{:?}",domain_spec[i]);
    }
 //   pretty_print_int(&mut domain_spec );
   
   // pretty_print_vec(&mut p);
    println!("-------------------\n");


 //   pretty_print_vec(&mut b);
    'running: for k in 0..100{
        let ten_millis = time::Duration::from_millis(100);
        //std::thread::sleep(ten_millis);
        canvas.set_draw_color(Color::RGB(0,0,0));
        canvas.present();
        canvas.clear();
      // canvas.set_draw_color(Color::RGB(0,0,0));

       // canvas.fill_rect(Rect::new(0, 0, 800, 600));
       for event in event_pump.poll_iter() {
        match event {
            Event::Quit { .. }
            | Event::KeyDown {
                keycode: Some(Keycode::Escape),
                ..
            } => break 'running,
            _ => {}
        }}
       for i in 0..u.len(){

            for j in 0..u[0].len(){

                if domain_spec[i][j]!=-1{

                    let grayscale = u[i][j] as u8;
                   //println!("{}", 2*(i as i32));
                    canvas.set_draw_color(Color::RGB(2*grayscale,2*grayscale, 2*grayscale));
                    canvas.fill_rect(Rect::new(2*(i as i32), 2*(j as i32), 2, 2));

                }

            }

        }
        
       // pretty_print_vec(&mut u);
/*
        ax = multiply_by_a_matrix(&mut u,&mut domain_spec,alpha_sim_par,gamma_sim_par);
        // pretty_print_vec(&mut ax);
        // println!("-------------------\n");
         b = compute_b(&mut u, &mut domain_spec, &mut control_array,gamma_sim_par, alpha_sim_par);
         r = subtrac_vec(&mut b, &mut domain_spec, &mut ax);
         p = force_copy(&mut r);
         delta_solve = scalar_product_itself(&mut r, &mut domain_spec );
*/
         if k == 50{
            control_array = vec![0.0; 2*A-2*B];
        }
    }

    Ok(())
}
