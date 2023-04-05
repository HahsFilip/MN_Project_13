use rand::Rng;
extern crate sdl2;
use colors_transform;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::rect::Rect;
//use colors_transform::Color;
//use colors_transform::{Hsl, Color};
use colors_transform::Color as OtherColor;
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
fn multiply_by_a_matrix_star( x: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>, alpha: f32, gamma: f32 ) -> Vec<Vec<f32>> {
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
fn compute_b_star (u: &mut Vec<Vec<f32>>,domain: &mut Vec<Vec<i32>>, c:&mut Vec<f32>, gamma: f32, alpha: f32 ) -> Vec<Vec<f32>>{
    let mut result = vec![vec![0.0; u[0].len()]; u.len()];
    for i in 0..u.len(){
        for j in 0..u[0].len(){
            if domain[i][j] != -1{
                if domain[i][j] == -2{
                    result[i][j] = u[i][j];
                
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
fn average_of_vec (x:  Vec<Vec<f32>>,domain:  Vec<Vec<i32>>)-> f32{
        let mut result: f32 = 0.0;
        let mut it = 0.0;
        for i in 0..x.len(){
            for j in 0..x[0].len(){
                if domain[i][j] != -1{
                    result = result + x[i][j];
                    it = it +1.0;
                }
            }
        }
        return result/it;
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

fn conjugate_gradiant(gamma_sim: f32,  alpha_sim: f32 ,
    u_0:  Vec<Vec<f32>>,domain_spec: &mut Vec<Vec<i32>>,
    a_func: fn( &mut Vec<Vec<f32>>, &mut Vec<Vec<i32>>, f32, f32) -> Vec<Vec<f32>>,
    b_func: fn( &mut Vec<Vec<f32>>, &mut Vec<Vec<i32>>,&mut Vec<f32>, f32,f32) -> Vec<Vec<f32>>,
    c_array:  Option<Vec<f32>> ) -> Vec<Vec<f32>>{
   
   
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
    ax = a_func( &mut u,domain_spec,alpha_sim,gamma_sim);
    // pretty_print_vec(&mut ax);
    // println!("-------------------\n");
     b = b_func(&mut u,  domain_spec, &mut c_array.unwrap_or(vec![1.1; 10]),gamma_sim, alpha_sim);
     r = subtrac_vec( &mut b, domain_spec, &mut ax);
     p = r.clone();
     delta_solve = scalar_product_itself(&mut r, domain_spec );
     let gamma_zero = delta_solve;
     for _n in 0..100{
        z = multiply_by_a_matrix(&mut p, domain_spec,alpha_sim, gamma_sim);
        
        alpha_solve = -delta_solve/scalar_product(&mut p,  domain_spec, &mut z);
       
        tmp = multiply_by_scalar_vec(&mut p,  domain_spec, alpha_solve);
       
        u = subtrac_vec(&mut u, domain_spec, &mut tmp);

        tmp = multiply_by_scalar_vec(&mut z, domain_spec,-alpha_solve);
        r = subtrac_vec(&mut r,  domain_spec, &mut tmp);
        gamma_solve = scalar_product_itself(&mut r, domain_spec);
       // println!("{}", gamma_solve/gamma_zero);
        let distance = scalar_product_itself(&mut r,  domain_spec);

        beta_solve = gamma_solve/delta_solve;
        tmp = multiply_by_scalar_vec(&mut p,  domain_spec, -beta_solve);
        p = subtrac_vec(&mut r, domain_spec, &mut tmp);
        delta_solve = gamma_solve;
        if gamma_solve/gamma_zero < 0.000001{
            break;
        }
       
    }
     u
}

fn change_control_array( adjoint_matrix:  Vec<Vec<f32>>,  domain:  Vec<Vec<i32>>, c_array: &mut Vec<f32>, alpha: f32, eps: f32) -> Vec<f32>{
    let mut result = vec![25.0;c_array.len()];
    for i in 0..adjoint_matrix.len(){
        for j in 0..adjoint_matrix[0].len(){
            if domain[i][j] != -1{
                if domain[i][j] != -2{   
                    let index: usize = domain[i][j] as usize;
                    result[index]= c_array[index]+ alpha*eps*adjoint_matrix[i][j];
                    }
            }
        
    }
}
result
}
fn main()-> Result<(), String> {
    let diffusivity = 10.0;
    let h = 1.0;
    let beta = 10.0;
    let dt = 0.010;
    let n_time_steps = 100;
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


    const A: usize = 40;
    const B: usize = 10;
    const C: usize = 30;
    const D: usize = 16;
    let a_int: i32 = A.try_into().unwrap();
    let b_int: i32 = B.try_into().unwrap();
    let c_int: i32 = C.try_into().unwrap();
    let d_int: i32 = D.try_into().unwrap();


    let mut domain_spec= vec![vec![-1; C+2]; A+2];
    let mut rng = rand::thread_rng();

    let mut u= vec![vec![1.1; C+2]; A+2];
    let mut u_backup = vec![vec![1.1; C+2]; A+2];
    let gamma_sim_par = 2.0*h*beta/diffusivity;
    let alpha_sim_par = diffusivity*dt/(h*h);
    let adjoint_gamma = 2.0*h*beta;

    for i in 0..A{
        for j in 0..C{
            domain_spec[i+1][j+1] = internal_detection(i.try_into().unwrap(), j.try_into().unwrap(), a_int,b_int,c_int,d_int);
            u[i+1][j+1] = rng.gen::<f32>()*50.1;
           //  if (i as i32 -20).abs() < 3{
           //      u[i+1][j+1]= 100.0;
          //   }
        }
       // println!("{:?}",domain_spec[i]);
    }
    let goal = average_of_vec(u.clone(), domain_spec.clone());
    let mut control_array = vec![vec![goal; 2*A-2*B];n_time_steps];

    u_backup = u.clone();
 //   pretty_print_int(&mut domain_spec );
   
   // pretty_print_vec(&mut p);
    println!("-------------------\n");

    for l in 0..500000-1{
        u = u_backup.clone();
 //   pretty_print_vec(&mut b);
    'running: for k in 0..n_time_steps{
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
            | Event::KeyDown     {
                keycode: Some(Keycode::Escape),
                ..
            } => break 'running,
            _ => {}
        }
    }
{
    let pixel_size: u32 = 10;
    let max_x: i32 = (pixel_size as i32)*(u.len() as i32);
    let max_y: i32 =(pixel_size as i32)*(u[0].len() as i32);
    let mut  max_temp_domain = goal;
    let mut  min_temp_domain = goal;
    for range_finder in 0.. u.len(){
        for range_finder_2 in 0.. u[range_finder].len(){
            if domain_spec[range_finder][range_finder_2] != -1{
                if max_temp_domain < u[range_finder][range_finder_2]{
                    max_temp_domain = u[range_finder][range_finder_2];
                }
                if min_temp_domain> u[range_finder][range_finder_2]{
                    min_temp_domain = u[range_finder][range_finder_2];
                }
            }
   

        }
    }
    let d_t_domain = max_temp_domain - min_temp_domain;
    let scale_domain = 350.0/d_t_domain;
    for i in 0..u.len(){

            for j in 0..u[0].len(){

                if domain_spec[i][j]!=-1{

                    let grayscale = ((u[i][j] -  min_temp_domain)*scale_domain) as u16;
                    let hex_color = colors_transform::Hsl::from((360-grayscale).into(), 100.0, 50.0);
                    let rgb = hex_color.to_rgb();
                    //println!("{}", 2*(i as i32))
                    canvas.set_draw_color(Color::RGB(rgb.get_red() as u8,rgb.get_green() as u8,rgb.get_blue() as u8));
                    canvas.fill_rect(Rect::new((pixel_size as i32)*(i as i32), (pixel_size as i32)*(j as i32)+pixel_size as i32, pixel_size , pixel_size) );
                    
                
                }

            }

        }
        let mut max_temp = control_array[0][0];
        let mut min_temp = control_array[0][0];
        for range_finder in 0.. control_array.len(){
            for range_finder_2 in 0.. control_array[range_finder].len(){
                if max_temp < control_array[range_finder][range_finder_2]{
                    max_temp = control_array[range_finder][range_finder_2];
                }
                if min_temp> control_array[range_finder][range_finder_2]{
                    min_temp = control_array[range_finder][range_finder_2];
                }

            }
        }
        let d_t = max_temp - min_temp;
        let scale = 200.0/d_t;
        for i in 0..control_array[k].len(){
           // println!("{}", d_t);
            let grayscale = ((control_array[k][i]-min_temp)*scale) as u8;
            //println!("{}", grayscale);
            canvas.set_draw_color(Color::RGB(grayscale,grayscale, grayscale));
   
            if i < (A-2*B){            
                canvas.fill_rect(Rect::new((pixel_size as i32)*((i +B +1)as i32),0 , pixel_size , pixel_size ));
            }else{
                canvas.fill_rect(Rect::new((pixel_size as i32)*(i as i32 - (A-2*B -1 )  as i32),max_y+pixel_size as i32, pixel_size , pixel_size ));

            }
            //println!("{}", 2*(i as i32));

             

        }
    }
        u = conjugate_gradiant(gamma_sim_par, alpha_sim_par, u, &mut domain_spec, multiply_by_a_matrix, compute_b, Some(control_array[k].clone()));
        

    }
    
    let mut u_star = vec![vec![goal; C+2]; A+2];
    u_star = subtrac_vec(&mut u_star,&mut domain_spec, &mut u);
    //pretty_print_vec(&mut u_star);
    
   // pretty_print_vec(&mut control_array);
    for k in 0..control_array.len(){
        let max_ind = control_array.len();
        u_star = conjugate_gradiant(adjoint_gamma, alpha_sim_par, u_star, &mut domain_spec, multiply_by_a_matrix_star, compute_b_star, None);
        control_array[n_time_steps-k-1] = change_control_array(u_star.clone(), domain_spec.clone(), &mut control_array[max_ind-1-k], beta, 1.1);
    }
    println!("{}", scalar_product_itself( &mut u_star, &mut  domain_spec));
}
    Ok(())
}
