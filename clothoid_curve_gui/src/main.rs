//! Adapted from egui custom_plot_manipulation example
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use clothoid_curve::Clothoid;
use eframe::egui::{self, DragValue, Event, Vec2};
use egui_plot::{Legend, Line, LineStyle, PlotPoints};


fn main() -> Result<(), eframe::Error> {
    env_logger::init(); // Log to stderr (if you run with `RUST_LOG=debug`).
    let options = eframe::NativeOptions::default();
    eframe::run_native(
        "CurvePlot",
        options,
        Box::new(|_cc| Box::<PlotCurve>::default()),
    )
}

struct PlotCurve {
    lock_x: bool,
    lock_y: bool,
    ctrl_to_zoom: bool,
    shift_to_horizontal: bool,
    zoom_speed: f32,
    scroll_speed: f32,
    x0: f64,
    y0: f64,
    theta0: f64,
    kappa0: f64,
    dk: f64,
    length: f64,
    count: i64,
}

impl Default for PlotCurve {
    fn default() -> Self {
        Self {
            lock_x: false,
            lock_y: false,
            ctrl_to_zoom: false,
            shift_to_horizontal: false,
            zoom_speed: 1.0,
            scroll_speed: 1.0,
            x0: 0.0,
            y0: 0.0,
            theta0: 0.0,
            kappa0: 0.0,
            dk: 5.0,
            length: 1.0,
            count: 0,
        }
    }
}

impl eframe::App for PlotCurve {
    fn update(&mut self, ctx: &egui::Context, _: &mut eframe::Frame) {
        egui::SidePanel::left("options").show(ctx, |ui| {

            ui.horizontal(|ui| {
                ui.add(
                    DragValue::new(&mut self.x0)
                        .clamp_range(-10.0..=10.0)
                        .speed(0.1),
                );
                ui.label("x0").on_hover_text("x0");
            });

            ui.horizontal(|ui| {
                ui.add(
                    DragValue::new(&mut self.y0)
                        .clamp_range(-10.0..=10.0)
                        .speed(0.1),
                );
                ui.label("y0").on_hover_text("y0");
            });

            ui.horizontal(|ui| {
                ui.add(
                    DragValue::new(&mut self.theta0)
                        .clamp_range(-10.0..=10.0)
                        .speed(0.1),
                );
                ui.label("theta0").on_hover_text("theta0");
            });

            ui.horizontal(|ui| {
                ui.add(
                    DragValue::new(&mut self.kappa0)
                        .clamp_range(-2.0..=2.0)
                        .speed(0.1),
                );
                ui.label("kappa0").on_hover_text("kappa0");
            });

            ui.horizontal(|ui| {
                ui.add(
                    DragValue::new(&mut self.dk)
                        .clamp_range(-10.0..=10.0)
                        .speed(0.01),
                );
                ui.label("dk").on_hover_text("dk");
            });

            ui.horizontal(|ui| {
                ui.add(
                    DragValue::new(&mut self.length)
                        .clamp_range(0.1..=20.0)
                        .speed(0.1),
                );
                ui.label("length").on_hover_text("length");
            });

            let x_text = "Check to keep the X axis fixed, i.e., pan and zoom will only affect the Y axis";
            ui.checkbox(&mut self.lock_x, "Lock x axis").on_hover_text(x_text);
            let y_text = "Check to keep the Y axis fixed, i.e., pan and zoom will only affect the X axis";
            ui.checkbox(&mut self.lock_y, "Lock y axis").on_hover_text(y_text);
            let ctrl_text = "If unchecked, the behavior of the Ctrl key is inverted compared to the default controls\ni.e., scrolling the mouse without pressing any keys zooms the plot";
            ui.checkbox(&mut self.ctrl_to_zoom, "Ctrl to zoom").on_hover_text(ctrl_text);
            let shift_text = "If unchecked, the behavior of the shift key is inverted compared to the default controls\ni.e., hold to scroll vertically, release to scroll horizontally";
            ui.checkbox(&mut self.shift_to_horizontal, "Shift for horizontal scroll").on_hover_text(shift_text);
            ui.horizontal(|ui| {
                ui.add(
                    DragValue::new(&mut self.zoom_speed)
                        .clamp_range(0.1..=2.0)
                        .speed(0.1),
                );
                ui.label("Zoom speed").on_hover_text("How fast to zoom in and out with the mouse wheel");
            });
            ui.horizontal(|ui| {
                ui.add(
                    DragValue::new(&mut self.scroll_speed)
                        .clamp_range(0.1..=100.0)
                        .speed(0.1),
                );
                ui.label("Scroll speed").on_hover_text("How fast to pan with the mouse wheel");
            });
        });
        egui::CentralPanel::default().show(ctx, |ui| {
            let (scroll, pointer_down, modifiers) = ui.input(|i| {
                let scroll = i.events.iter().find_map(|e| match e {
                    Event::MouseWheel {
                        unit: _,
                        delta,
                        modifiers: _,
                    } => Some(*delta),
                    _ => None,
                });
                (scroll, i.pointer.primary_down(), i.modifiers)
            });

            ui.label("This example shows how to use raw input events to implement different plot controls than the ones egui provides by default, e.g., default to zooming instead of panning when the Ctrl key is not pressed, or controlling much it zooms with each mouse wheel step.");

            egui_plot::Plot::new("plot")
                .allow_zoom(false)
                .allow_drag(false)
                .allow_scroll(false)
                .legend(Legend::default())
                .show(ui, |plot_ui| {
                    if let Some(mut scroll) = scroll {
                        if modifiers.ctrl == self.ctrl_to_zoom {
                            scroll = Vec2::splat(scroll.x + scroll.y);
                            let mut zoom_factor = Vec2::from([
                                (scroll.x * self.zoom_speed / 10.0).exp(),
                                (scroll.y * self.zoom_speed / 10.0).exp(),
                            ]);
                            if self.lock_x {
                                zoom_factor.x = 1.0;
                            }
                            if self.lock_y {
                                zoom_factor.y = 1.0;
                            }
                            plot_ui.zoom_bounds_around_hovered(zoom_factor);
                        } else {
                            if modifiers.shift == self.shift_to_horizontal {
                                scroll = Vec2::new(scroll.y, scroll.x);
                            }
                            if self.lock_x {
                                scroll.x = 0.0;
                            }
                            if self.lock_y {
                                scroll.y = 0.0;
                            }
                            let delta_pos = self.scroll_speed * scroll;
                            plot_ui.translate_bounds(delta_pos);
                        }
                    }
                    if plot_ui.response().hovered() && pointer_down {
                        let mut pointer_translate = -plot_ui.pointer_coordinate_drag_delta();
                        if self.lock_x {
                            pointer_translate.x = 0.0;
                        }
                        if self.lock_y {
                            pointer_translate.y = 0.0;
                        }
                        plot_ui.translate_bounds(pointer_translate);
                    }

                    let curve = Clothoid::create(
                        self.x0,
                        self.y0,
                        self.theta0,
                        self.kappa0,
                        self.dk,
                        self.length,
                    );
                    let xy = curve.get_points(100);
                    /*
                    if self.count % 100 == 0 {
                        println!("{}",  self.dk);
                        println!("{:?}", xy);
                    }
                    */
                    self.count += 1;
                    let curve_points = PlotPoints::new(xy);
                    plot_ui.line(Line::new(curve_points).name("Curve").style(LineStyle::dashed_dense()));
                });
        });
    }
}
