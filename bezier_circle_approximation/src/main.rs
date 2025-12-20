/*!
Lucas Walter

December 2025

show how well a set of bezier curves can approximate a circle
*/
use egui::{CentralPanel, Color32, Stroke, TopBottomPanel};
use egui_plot::{Legend, Line};
use std::f64::consts::PI;
use stroke::f64::{CubicBezier, Point, PointN};
// use tracing::{debug, error, info, warn};

struct BezierCircleApproximation {
    radius: f64,
    num: usize,
    angle: f64,
    handle_length: f64,
}

impl eframe::App for BezierCircleApproximation {
    fn update(&mut self, ctx: &egui::Context, _: &mut eframe::Frame) {
        let mut bezier_distance = Vec::new();
        TopBottomPanel::top("controls").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("radius");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.radius)
                        .speed(0.003)
                        .range(0.2..=4.0)
                        .update_while_editing(false),
                );

                ui.label("bezier handle length");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.handle_length)
                        .speed(0.003)
                        .range(0.01..=5.0)
                        .update_while_editing(false),
                );

                // https://stackoverflow.com/a/27863181/603653
                let num = 2.0 * PI / self.angle;
                let optimal_length = self.radius * 4.0 / 3.0 * (PI / (2.0 * num)).tan();

                if ui
                    .button(format!("optimal {:.6}", optimal_length))
                    .clicked()
                {
                    self.handle_length = optimal_length;
                }

                ui.label("num segments");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.num)
                        .speed(0.04)
                        .range(1..=16)
                        .update_while_editing(false),
                );

                ui.label("segment angle");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.angle)
                        .speed(0.01)
                        .range(0.1..=2.0 * PI)
                        .update_while_editing(false),
                );
            });
        });

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
                    .show(ui, |plot_ui| {
                        let mut start_angle: f64 = 0.0;
                        for i in 0..self.num {
                            let mut circle_pts = Vec::new();
                            let max_angle_ind = 32;
                            let fr = self.angle / max_angle_ind as f64;

                            let mut angle = start_angle;

                            // the bezier end points
                            let bz_pt0: [f64; 2] =
                                [self.radius * angle.cos(), self.radius * angle.sin()];
                            // angles from the end points to the handles
                            let bz_angle0 = angle + PI / 2.0;

                            for angle_ind in 0..(max_angle_ind + 1) {
                                let x = self.radius * angle.cos();
                                let y = self.radius * angle.sin();
                                circle_pts.push([x, y]);
                                if angle_ind < max_angle_ind {
                                    angle += fr;
                                }
                            }

                            let bz_pt3 = [self.radius * angle.cos(), self.radius * angle.sin()];
                            let bz_angle1 = angle - PI / 2.0;

                            let color = {
                                // if i.is_multiple_of(2) {
                                if i % 2 == 0 {
                                    Color32::PURPLE
                                } else {
                                    Color32::BLUE
                                }
                            };
                            plot_ui.line(
                                Line::new(circle_pts)
                                    .name("segment".to_string())
                                    .allow_hover(false)
                                    .stroke(Stroke::new(2.0, color)),
                            );

                            // find where the handles are
                            let bz_pt1 = [
                                bz_pt0[0] + self.handle_length * bz_angle0.cos(),
                                bz_pt0[1] + self.handle_length * bz_angle0.sin(),
                            ];
                            let bz_pt2 = [
                                bz_pt3[0] + self.handle_length * bz_angle1.cos(),
                                bz_pt3[1] + self.handle_length * bz_angle1.sin(),
                            ];
                            let bezier: CubicBezier<PointN<2>, 2> = CubicBezier::new(
                                PointN::new(bz_pt0),
                                PointN::new(bz_pt1),
                                PointN::new(bz_pt2),
                                PointN::new(bz_pt3),
                            );

                            plot_ui.line(
                                Line::new(vec![bz_pt0, bz_pt1])
                                    .name("handle".to_string())
                                    .allow_hover(false)
                                    .stroke(Stroke::new(2.0, Color32::CYAN)),
                            );

                            plot_ui.line(
                                Line::new(vec![bz_pt3, bz_pt2])
                                    .name("handle".to_string())
                                    .allow_hover(false)
                                    .stroke(Stroke::new(2.0, Color32::MAGENTA)),
                            );

                            let bezier_length = bezier.arclen_castlejau(None);

                            let mut bezier_pts = Vec::new();
                            let num_t = 32;
                            let fr = 1.0 / num_t as f64;
                            let mut euclidean_tfrac = 0.0;
                            for _i in 0..(num_t + 1) {
                                let desired_length = euclidean_tfrac * bezier_length;
                                let (_len, parametric_tfrac) =
                                    bezier.desired_len_to_parametric_t(desired_length, None);
                                let pt = bezier.eval(parametric_tfrac);
                                let pt = [pt.axis(0), pt.axis(1)];
                                // TODO(lucasw) get
                                let dist = (pt[0] * pt[0] + pt[1] * pt[1]).sqrt();
                                let angle = pt[1].atan2(pt[0]);
                                // bezier_distance.push([angle, dist, self.radius]);
                                bezier_distance.push([angle, dist, self.radius]);
                                // info!("{angle:.3}, {dist:.3}");
                                // how far it is from 1.0 is the bezier approximation error
                                bezier_pts.push(pt);
                                euclidean_tfrac += fr;
                            }

                            let color = {
                                // if i.is_multiple_of(2) {
                                if i % 2 == 0 {
                                    Color32::GREEN
                                } else {
                                    Color32::LIGHT_GREEN
                                }
                            };

                            plot_ui.points(
                                egui_plot::Points::new(bezier_pts.clone())
                                    .name("bezier points")
                                    .allow_hover(false)
                                    .radius(2.5)
                                    .color(Color32::GOLD),
                            );

                            plot_ui.line(
                                Line::new(bezier_pts)
                                    .name("bezier".to_string())
                                    .allow_hover(false)
                                    .stroke(Stroke::new(2.0, color)),
                            );

                            start_angle += self.angle;
                        }
                    });

                ui.end_row();

                let distance_error: Vec<[f64; 2]> = bezier_distance
                    .clone()
                    .into_iter()
                    .map(|x| [x[0], x[1] - x[2]])
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
                            egui_plot::Points::new(distance_error)
                                .name("bezier points")
                                .allow_hover(false)
                                .radius(2.5)
                                .color(Color32::GOLD),
                        );
                    });
            });
        });

        TopBottomPanel::bottom("stats").show(ctx, |ui| {
            let mut measured = Vec::new();
            let mut expected = Vec::new();
            for x in bezier_distance {
                measured.push(x[1]);
                expected.push(x[2]);
            }
            // println!("{vec_vec:?}");
            let rmse = eval_metrics::regression::rmse(&measured, &expected).unwrap();
            // ui.label(format!("rmse {rmse:.6}"));
            ui.label(format!("rmse {rmse:.9}"));
        });
    }
}

impl BezierCircleApproximation {
    fn new(_cc: &eframe::CreationContext<'_>) -> Result<Self, anyhow::Error> {
        Ok(BezierCircleApproximation {
            radius: 1.0,
            num: 4,
            angle: PI / 2.0,
            handle_length: 0.4,
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
        Box::new(|cc| Ok(Box::new(BezierCircleApproximation::new(cc)?))),
    );

    Ok(())
}
