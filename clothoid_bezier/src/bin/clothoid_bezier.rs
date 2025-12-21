/*!
Lucas Walter

December 2025

show how well a set of bezier curves can approximate a clothoid curve
*/
use clothoid_curve::clothoid::Clothoid;
use egui::{CentralPanel, Color32, Stroke, TopBottomPanel};
use egui_plot::{Legend, Line, Points};
use std::f64::consts::PI;
use stroke::f64::{CubicBezier, Point, PointN};
// use tracing::{debug, error, info, warn};

struct ClothoidBezierApproximation {
    curvature: f64,
    curvature_rate: f64,
    length: f64,
    // num: usize,
    handle_length0: f64,
    handle_length1: f64,
    // going_to_optimal: bool,
}

/*
// TODO(lucasw) put in lib.rs, add unit tests
/// https://stackoverflow.com/a/27863181/603653
/// this will work poorly for small curvatures
fn radius_to_bezier_handle_length(radius: f64, arc_angle: f64) -> f64 {
    // TODO(lucasw) check if handle length is > arc_len / 2.0 or similar
    // let arc_len = radius * arc_angle;
    // let num = 2.0 * PI / arc_angle;
    // radius * 4.0 / 3.0 * (PI / (2.0 * num)).tan()
    radius * 4.0 / 3.0 * (arc_angle / 4.0).tan()
}
*/

/*
fn curvature_to_bezier_handle_length(curvature: f64, arc_angle: f64) -> f64 {
    radius_to_bezier_handle_length(1.0 / curvature, arc_angle)
}
*/

impl eframe::App for ClothoidBezierApproximation {
    fn update(&mut self, ctx: &egui::Context, _: &mut eframe::Frame) {
        // let mut bezier_distance = Vec::new();
        TopBottomPanel::top("controls").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("curvature");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.curvature)
                        .speed(0.003)
                        .range(-4.0..=4.0)
                        .update_while_editing(false),
                );

                ui.label("curvature rate");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.curvature_rate)
                        .speed(0.003)
                        .range(-4.0..=4.0)
                        .update_while_editing(false),
                );

                ui.label("length");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.length)
                        .speed(0.01)
                        .range(0.1..=20.0)
                        .update_while_editing(false),
                );

                ui.label("bezier handle length");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.handle_length0)
                        .speed(0.003)
                        .range(0.0..=5.0)
                        .update_while_editing(false),
                );

                let _resp = ui.add(
                    egui::DragValue::new(&mut self.handle_length1)
                        .speed(0.003)
                        .range(0.0..=5.0)
                        .update_while_editing(false),
                );

                /*
                let optimal_length = radius_to_bezier_handle_length(self.radius, self.angle);

                if ui
                    .button(format!("optimal {:.6}", optimal_length))
                    // .is_pointer_button_down_on()
                    .clicked()
                {
                    // go to optimal in one step
                    // self.handle_length = optimal_length;
                    self.going_to_optimal = true;
                }

                if self.going_to_optimal {
                    // smoothly move towards optimal
                    // info!("pressed");
                    let error = (optimal_length - self.handle_length).clamp(-0.015, 0.015);
                    self.handle_length += error;
                    // if nothing changes won't get a repaint so won't get to optimal,
                    // so force a repaint so this keeps executing
                    if error != 0.0 {
                        ctx.request_repaint();
                    } else {
                        self.going_to_optimal = false;
                    }
                }
                */
            });
            /*
            ui.horizontal(|ui| {
                let full_circle_angle = 2.0 * PI / (self.num as f64);
                if ui
                    .button(format!("full circle {:.2}", full_circle_angle))
                    .clicked()
                {
                    self.angle = full_circle_angle;
                }
            });
            */
        });

        let mut bezier_distance = Vec::new();
        CentralPanel::default().show(ctx, |ui| {
            egui::Grid::new("grid").num_columns(1).show(ui, |ui| {
                egui_plot::Plot::new("plot")
                    .auto_bounds(true)
                    .allow_double_click_reset(true)
                    .allow_zoom(true)
                    .allow_drag(true)
                    .allow_scroll(true)
                    .legend(Legend::default())
                    .data_aspect(1.0)
                    .view_aspect(1.0)
                    // .show_grid(false)
                    .show(ui, |plot_ui| {
                        let x0 = 0.0;
                        let y0 = 0.0;

                        let num_pts = ((self.length * 10.0) as u32).clamp(4, 128);

                        let clothoid = Clothoid::create(
                            x0,
                            y0,
                            0.0,
                            self.curvature,
                            self.curvature_rate,
                            self.length,
                        );

                        let clothoid_pts = clothoid.get_points(num_pts);

                        plot_ui.line(
                            Line::new(clothoid_pts.clone())
                                .name("clothoid")
                                .allow_hover(false)
                                .stroke(Stroke::new(2.0, Color32::MAGENTA)),
                        );

                        plot_ui.points(
                            Points::new(clothoid_pts)
                                .name("clothoid pts")
                                .allow_hover(true)
                                .radius(1.0)
                                .color(Color32::CYAN),
                        );

                        let clothoid_end = clothoid.get_end_clothoid();
                        // the bezier end points
                        let bz_pt0: [f64; 2] = [clothoid.x0, clothoid.y0];
                        let bz_pt3: [f64; 2] = [clothoid_end.x0, clothoid_end.y0];

                        // find where the handles are
                        let bz_angle0 = clothoid.theta0;
                        let bz_pt1 = [
                            bz_pt0[0] + self.handle_length0 * bz_angle0.cos(),
                            bz_pt0[1] + self.handle_length0 * bz_angle0.sin(),
                        ];

                        let bz_angle1 = clothoid_end.theta0;
                        let bz_pt2 = [
                            bz_pt3[0] - self.handle_length1 * bz_angle1.cos(),
                            bz_pt3[1] - self.handle_length1 * bz_angle1.sin(),
                        ];

                        let bezier: CubicBezier<PointN<2>, 2> = CubicBezier::new(
                            PointN::new(bz_pt0),
                            PointN::new(bz_pt1),
                            PointN::new(bz_pt2),
                            PointN::new(bz_pt3),
                        );

                        plot_ui.line(
                            Line::new(vec![bz_pt0, bz_pt1])
                                .name("handle")
                                .allow_hover(false)
                                .stroke(Stroke::new(2.0, Color32::CYAN)),
                        );

                        plot_ui.line(
                            Line::new(vec![bz_pt3, bz_pt2])
                                .name("handle")
                                .allow_hover(false)
                                .stroke(Stroke::new(2.0, Color32::MAGENTA)),
                        );

                        plot_ui.points(
                            Points::new(vec![bz_pt1, bz_pt2])
                                .name("handle")
                                .allow_hover(true)
                                .radius(4.0)
                                .color(Color32::CYAN),
                        );

                        // TODO(lucasw) add this to stroke 'get_euclidean_pts(num)'
                        let bezier_length = bezier.arclen_castlejau(None);
                        let fr = 1.0 / (num_pts - 1) as f64;

                        let mut bezier_pts = Vec::new();
                        for i in 0..num_pts {
                            let euclidean_tfrac = i as f64 * fr;
                            let s_distance = euclidean_tfrac * bezier_length;

                            let desired_length = euclidean_tfrac * bezier_length;
                            let (_len, parametric_tfrac) =
                                bezier.desired_len_to_parametric_t(desired_length, None);
                            let pt = bezier.eval(parametric_tfrac);
                            let pt = [pt.axis(0), pt.axis(1)];
                            bezier_pts.push(pt);

                            // this will be increasingly different from bezier s_distance
                            // the more the two curves differ
                            let clothoid_s_distance = euclidean_tfrac * clothoid.length;
                            let cpt = clothoid.get_clothoid(clothoid_s_distance);
                            let ex = pt[0] - cpt.x0;
                            let ey = pt[1] - cpt.y0;
                            let error_distance = (ex * ex + ey * ey).sqrt();

                            let curvature = bezier.curvature(parametric_tfrac);

                            bezier_distance.push((
                                s_distance,
                                error_distance,
                                curvature,
                                cpt.curvature(),
                            ));
                        }

                        let color = Color32::ORANGE;
                        plot_ui.points(
                            Points::new(bezier_pts.clone())
                                .name("bezier points")
                                .allow_hover(false)
                                .radius(2.0)
                                .color(color),
                        );

                        plot_ui.line(
                            Line::new(bezier_pts)
                                .name("bezier".to_string())
                                .allow_hover(false)
                                .stroke(Stroke::new(1.5, color)),
                        );
                    });

                ui.end_row();

                let position_error: Vec<[f64; 2]> = bezier_distance
                    .clone()
                    .into_iter()
                    .map(|x| [x.0, x.1])
                    .collect();
                let curvature_error: Vec<[f64; 2]> = bezier_distance
                    .clone()
                    .into_iter()
                    .map(|x| [x.0, x.2 - x.3])
                    .collect();

                egui_plot::Plot::new("error")
                    .auto_bounds(true)
                    .allow_double_click_reset(true)
                    .allow_zoom(true)
                    .allow_drag(true)
                    .allow_scroll(true)
                    .legend(Legend::default())
                    // .data_aspect(10.0)
                    .view_aspect(3.0)
                    .show(ui, |plot_ui| {
                        plot_ui.points(
                            Points::new(position_error.clone())
                                .name("postion error points")
                                .allow_hover(false)
                                .radius(2.5)
                                .color(Color32::GOLD),
                        );

                        plot_ui.line(
                            Line::new(position_error)
                                .name("position error")
                                .allow_hover(false)
                                .stroke(Stroke::new(1.2, Color32::GOLD)),
                        );

                        plot_ui.line(
                            Line::new(curvature_error)
                                .name("curvature error")
                                .allow_hover(false)
                                .stroke(Stroke::new(1.2, Color32::PURPLE)),
                        );
                    });
            });
        });

        TopBottomPanel::bottom("stats").show(ctx, |ui| {
            let mut measured = Vec::new();
            let mut expected = Vec::new();
            for vals in &bezier_distance {
                measured.push(0.0);
                expected.push(vals.1);
            }
            let rmse_position = eval_metrics::regression::rmse(&measured, &expected).unwrap();
            ui.label(format!("rmse position {rmse_position:.9}"));

            let mut measured = Vec::new();
            let mut expected = Vec::new();
            for vals in &bezier_distance {
                measured.push(0.0);
                expected.push(vals.2 - vals.3);
            }
            let rmse_curvature = eval_metrics::regression::rmse(&measured, &expected).unwrap();
            // ui.label(format!("rmse {rmse:.6}"));
            ui.label(format!("rmse curvature {rmse_curvature:.9}"));
        });
    }
}

impl ClothoidBezierApproximation {
    fn new(_cc: &eframe::CreationContext<'_>) -> Result<Self, anyhow::Error> {
        Ok(ClothoidBezierApproximation {
            curvature: 0.0,
            curvature_rate: 1.275,
            length: PI / 2.0,
            handle_length0: 0.4,
            handle_length1: 0.4,
            // going_to_optimal: false,
        })
    }
}

fn main() -> Result<(), anyhow::Error> {
    let subscriber = tracing_subscriber::fmt()
        .compact()
        .with_file(true)
        .with_line_number(true)
        .with_thread_ids(true)
        .with_target(true)
        .finish();
    tracing::subscriber::set_global_default(subscriber)?;

    let options = eframe::NativeOptions {
        viewport: egui::viewport::ViewportBuilder::default()
            .with_inner_size(egui::vec2(720.0, 1024.0)),
        ..Default::default()
    };
    let _ = eframe::run_native(
        "Bezier Circle Approximation",
        options,
        Box::new(|cc| Ok(Box::new(ClothoidBezierApproximation::new(cc)?))),
    );

    Ok(())
}
