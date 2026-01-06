/*!
Lucas Walter

December 2025

show how well a set of bezier curves can approximate a clothoid curve
*/
use clothoid_bezier::f64::{
    cubic_bezier::EuclideanTFrac, BezierToClothoid, ClothoidBezierApproximation, CubicBezier2,
};
use clothoid_curve::f64::{
    curvature_per_meter, curvature_per_meter_float, AngleCosSin, Clothoid, CurvaturePerLength,
    Position,
};
use egui::{CentralPanel, Color32, Stroke, TopBottomPanel};
use egui_plot::{Legend, Line, Points};
use std::f64::consts::PI;
// use tracing::{debug, error, info, warn};
use tracing::warn;

use uom::num_traits::Zero;
use uom::si::{
    angle::degree,
    curvature::radian_per_meter,
    f64::{Angle, Curvature, Length},
    length::meter,
};

#[derive(Default)]
enum Mode {
    Bezier,
    #[default]
    Clothoid,
}

struct ClothoidToBezier {
    mode: Mode,
    clothoid_bezier: ClothoidBezierApproximation,
    bezier: CubicBezier2,
    handle_length0: f64,
    handle_length1: f64,
    end_x: f64,
    end_y: f64,
    end_angle: f64,
    // for animation
    target_handle_length0: f64,
    target_handle_length1: f64,
    going_to_target: bool,
}

impl eframe::App for ClothoidToBezier {
    fn update(&mut self, ctx: &egui::Context, _: &mut eframe::Frame) {
        let mut find_clothoid = false;
        let mut find_handles = false;
        let clothoid = self.clothoid_bezier.to_clothoid();
        let num_pts = 2 * ((self.clothoid_bezier.length * 5.0) as usize).clamp(4, 64);
        let fr = 1.0 / (num_pts - 1) as f64;

        TopBottomPanel::top("controls").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("curvature");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.clothoid_bezier.curvature)
                        .speed(0.003)
                        .range(-4.0..=4.0)
                        .update_while_editing(false),
                );

                ui.label("curvature rate");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.clothoid_bezier.curvature_rate)
                        .speed(0.003)
                        .range(-4.0..=4.0)
                        .update_while_editing(false),
                );

                ui.label("length");
                let _resp = ui.add(
                    egui::DragValue::new(&mut self.clothoid_bezier.length)
                        .speed(0.01)
                        .range(0.1..=20.0)
                        .update_while_editing(false),
                );

                if ui.button("find clothoid").clicked() {
                    self.mode = Mode::Bezier;
                    // keep initial curvature fixed, but solve for clothoid curvature rate
                    // and length that make the clothoid end point close to the bezier end
                    // point
                    find_clothoid = true;
                }
            });

            ui.horizontal(|ui| {
                ui.label("bezier handle length");
                let resp = ui.add(
                    egui::DragValue::new(&mut self.handle_length0)
                        .speed(0.003)
                        .range(0.0..=5.0)
                        .update_while_editing(false),
                );
                if resp.clicked() {
                    self.going_to_target = false;
                }

                let resp = ui.add(
                    egui::DragValue::new(&mut self.handle_length1)
                        .speed(0.003)
                        .range(0.0..=5.0)
                        .update_while_editing(false),
                );
                if resp.clicked() {
                    self.going_to_target = false;
                }

                if ui.button("find handles").clicked() {
                    self.mode = Mode::Clothoid;
                    find_handles = true;
                }
            });

            ui.horizontal(|ui| {
                let resp = ui.add(
                    egui::DragValue::new(&mut self.end_x)
                        .speed(0.002)
                        .update_while_editing(false),
                );
                if resp.dragged() {
                    self.mode = Mode::Bezier;
                }

                let resp = ui.add(
                    egui::DragValue::new(&mut self.end_y)
                        .speed(0.002)
                        .update_while_editing(false),
                );
                if resp.dragged() {
                    self.mode = Mode::Bezier;
                }

                let resp = ui.add(
                    egui::DragValue::new(&mut self.end_angle)
                        .speed(0.3)
                        .update_while_editing(false),
                );
                if resp.dragged() {
                    self.mode = Mode::Bezier;
                }
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

        let handle0 = Length::new::<meter>(self.handle_length0);
        let handle1 = Length::new::<meter>(self.handle_length1);
        match &self.mode {
            Mode::Bezier => {
                let start = clothoid.xy0;
                let end = Position {
                    x: Length::new::<meter>(self.end_x),
                    y: Length::new::<meter>(self.end_y),
                };
                let clothoid_end = Clothoid::create(
                    end,
                    AngleCosSin::from_angle(Angle::new::<degree>(self.end_angle)),
                    Curvature::zero(),
                    CurvaturePerLength::zero(),
                    handle1,
                );
                let clothoid_start_straight = clothoid.zero_curvature();
                let ctrl1 = clothoid_start_straight.get_xy(handle0);
                let ctrl2 = clothoid_end.get_xy(-handle1);

                self.bezier = CubicBezier2::new(&start, &ctrl1, &ctrl2, &end);

                if find_clothoid {
                    let length = self.bezier.arclen_castlejau();
                    let curvature_per_length =
                        curvature_per_meter(self.clothoid_bezier.curvature_rate);
                    let mut bezier_to_clothoid = BezierToClothoid::from_bezier(self.bezier.clone());
                    bezier_to_clothoid
                        .start_clothoid
                        .set_curvature(Curvature::new::<radian_per_meter>(
                            self.clothoid_bezier.curvature,
                        ));

                    let rv = bezier_to_clothoid.find_clothoid(curvature_per_length, length);

                    if let Ok(clothoid) = rv {
                        self.clothoid_bezier = ClothoidBezierApproximation {
                            curvature: clothoid.curvature().get::<radian_per_meter>(),
                            curvature_rate: curvature_per_meter_float(clothoid.curvature_rate()),
                            length: clothoid.length.get::<meter>(),
                        };
                    }
                }
            }
            Mode::Clothoid => {
                self.bezier = ClothoidBezierApproximation::get_bezier(&clothoid, handle0, handle1);

                if find_handles {
                    let rv = self.clothoid_bezier.find_handles(handle0, handle1);
                    if let Ok(handles) = rv {
                        // go to optimal in one step
                        // self.handle_length0 = handles.0;
                        // self.handle_length1 = handles.1;
                        self.target_handle_length0 = handles.0;
                        self.target_handle_length1 = handles.1;
                        self.going_to_target = true;
                    } else {
                        warn!("{rv:?}");
                    }
                }

                if self.going_to_target {
                    // smoothly move towards optimal
                    // info!("pressed");
                    let error0 =
                        (self.target_handle_length0 - self.handle_length0).clamp(-0.015, 0.015);
                    self.handle_length0 += error0;
                    let error1 =
                        (self.target_handle_length1 - self.handle_length1).clamp(-0.015, 0.015);
                    self.handle_length1 += error1;
                    // if nothing changes won't get a repaint so won't get to optimal,
                    // so force a repaint so this keeps executing
                    if error0 != 0.0 || error1 != 0.0 {
                        ctx.request_repaint();
                    } else {
                        self.going_to_target = false;
                    }
                }

                let end_clothoid = clothoid.get_end_clothoid();
                let end_pos = end_clothoid.xy0;
                self.end_x = end_pos.x.get::<meter>();
                self.end_y = end_pos.y.get::<meter>();
                self.end_angle = end_clothoid.theta0.get::<degree>();
            }
        }

        #[derive(Clone)]
        struct ClothoidBezierDistance {
            bezier_s_distance: Length,
            error_distance: Length,
            error_curvature: Curvature,
            // bezier_curvature: Curvature,
            // clothoid_curvature: Curvature,
        }
        let mut bezier_distance: Vec<ClothoidBezierDistance> = Vec::new();
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
                        let clothoid_pts = clothoid.get_points_num(num_pts);

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

                        let bz_pt0 = self.bezier.start().as_array_meter();
                        let bz_pt1 = self.bezier.ctrl1().as_array_meter();
                        let bz_pt2 = self.bezier.ctrl2().as_array_meter();
                        let bz_pt3 = self.bezier.end().as_array_meter();
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
                        let bezier_length = self.bezier.arclen_castlejau();

                        let mut bezier_pts = Vec::new();
                        for i in 0..num_pts {
                            let desired_euclidean_tfrac = EuclideanTFrac(i as f64 * fr);
                            let (achieved_euclidean_tfrac, parametric_tfrac) = self
                                .bezier
                                .euclidean_to_parametric_t(desired_euclidean_tfrac, bezier_length);
                            let bezier_s_distance = achieved_euclidean_tfrac.0 * bezier_length;
                            let pt = self.bezier.eval(parametric_tfrac);
                            bezier_pts.push(pt.as_array_meter());

                            // this will be increasingly different from bezier s_distance
                            // the more the two curves differ in length
                            let clothoid_s_distance = achieved_euclidean_tfrac.0 * clothoid.length;
                            let cpt = clothoid.get_clothoid(clothoid_s_distance);
                            let ex = pt.x - cpt.xy0.x;
                            let ey = pt.y - cpt.xy0.y;
                            let error_distance = (ex * ex + ey * ey).sqrt();

                            let bezier_curvature = self.bezier.curvature(parametric_tfrac);
                            let clothoid_curvature = cpt.curvature();
                            let error_curvature: Curvature = bezier_curvature - clothoid_curvature; //.into();

                            bezier_distance.push(ClothoidBezierDistance {
                                bezier_s_distance,
                                error_distance,
                                error_curvature,
                            });
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
                    .map(|v| {
                        [
                            v.bezier_s_distance.get::<meter>(),
                            v.error_distance.get::<meter>(),
                        ]
                    })
                    .collect();
                let curvature_error: Vec<[f64; 2]> = bezier_distance
                    .clone()
                    .into_iter()
                    .map(|v| {
                        [
                            v.bezier_s_distance.get::<meter>(),
                            v.error_curvature.get::<radian_per_meter>(),
                        ]
                    })
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
                expected.push(vals.error_distance.get::<meter>());
            }
            let rmse_position = eval_metrics::regression::rmse(&measured, &expected).unwrap();
            ui.label(format!("rmse position {rmse_position:.9}"));

            let mut measured = Vec::new();
            let mut expected = Vec::new();
            for vals in &bezier_distance {
                measured.push(0.0);
                expected.push(vals.error_curvature.get::<radian_per_meter>());
            }
            let rmse_curvature = eval_metrics::regression::rmse(&measured, &expected).unwrap();
            // ui.label(format!("rmse {rmse:.6}"));
            ui.label(format!("rmse curvature {rmse_curvature:.9}"));
        });
    }
}

impl ClothoidToBezier {
    fn new(_cc: &eframe::CreationContext<'_>) -> Result<Self, anyhow::Error> {
        let clothoid_bezier = ClothoidBezierApproximation {
            curvature: 1.0,
            curvature_rate: 0.0, // 1.275,
            length: PI / 2.0,
        };
        let clothoid = clothoid_bezier.to_clothoid();
        let end_clothoid = clothoid.get_end_clothoid();
        let end_pos = end_clothoid.xy0;
        let bezier =
            ClothoidBezierApproximation::get_bezier(&clothoid, Length::zero(), Length::zero());

        Ok(Self {
            mode: Mode::Clothoid,
            clothoid_bezier,
            bezier,
            handle_length0: 0.0,
            handle_length1: 0.0,
            target_handle_length0: 0.4,
            target_handle_length1: 0.4,
            going_to_target: true,
            end_x: end_pos.x.get::<meter>(),
            end_y: end_pos.y.get::<meter>(),
            end_angle: end_clothoid.theta0.get::<degree>(),
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
        "Clothoid Bezier Approximation",
        options,
        Box::new(|cc| Ok(Box::new(ClothoidToBezier::new(cc)?))),
    );

    Ok(())
}
