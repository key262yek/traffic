use moldybrody::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, State)]
pub struct Car {
    pos: Cartessian1D<f64>,
    vel: Cartessian1D<f64>,
    lane: usize,
    size: f64,
    max_speed: f64,
    drift: f64,
    mass: f64,
    behavior: f64,
}

impl Car {
    pub fn new(lane: usize, size: f64, max_speed: f64, drift: f64, behavior: f64) -> Self {
        Self {
            pos: Cartessian1D { coord: [0f64] },
            vel: Cartessian1D { coord: [0f64] },
            lane,
            size,
            max_speed,
            drift,
            mass: 1f64,
            behavior,
        }
    }

    pub fn drift_force(&self) -> Cartessian1D<f64> {
        Cartessian1D::new([self.drift * (self.max_speed - self.vel[0]).tanh()])
    }

    pub fn set_max_speed(&mut self, speed_limit: f64) {
        self.max_speed = (1.0 + self.behavior) * speed_limit;
    }

    pub fn safe_distance(&self, dir: bool) -> f64 {
        // dir == true : in front of self
        let v = self.vel[0];
        if dir {
            self.size + v * (2.0 + 0.05 * v)
        } else {
            self.size
        }
    }

    pub fn force_from(&self, other: &Self) -> Cartessian1D<f64> {
        // dir == true : other is in front of self
        let r = (&self.pos).distance(&other.pos);
        if self.pos[0] < other.pos[0] {
            let sd = self.safe_distance(true);
            let t = sd / r;
            Cartessian1D::new([-4.0 / r * (12.0 * t.powi(12) - 6.0 * t.powi(6))])
        } else {
            let sd = self.safe_distance(false);
            let t = sd / r;
            Cartessian1D::new([4.0 / r * (12.0 * t.powi(12) - 6.0 * t.powi(6))])
        }
    }
}

pub struct SpeedCam {
    pos: Cartessian1D<f64>,
    speed_limit: f64,
}

impl SpeedCam {
    pub fn new(pos: Cartessian1D<f64>, speed_limit: f64) -> Self {
        Self { pos, speed_limit }
    }

    pub fn set_max_speed(&self, car: &mut Car) {
        car.set_max_speed(self.speed_limit);
    }

    pub fn force_to(&self, car: &Car) -> Cartessian1D<f64> {
        if car.vel[0] < self.speed_limit {
            Cartessian1D::zeros()
        } else {
            Cartessian1D::zeros()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_car_acceleration() {
        let mut timeiter = ConstStep::<f64>::new(1e-1).unwrap();
        timeiter.set_tmax(20f64).unwrap();
        let mut car = Car::new(0, 1f64, 10f64, 1f64, 0f64);

        for (_t, dt) in timeiter.into_diff() {
            let force = car.drift_force();
            let movement = car.euler(&force, dt);
            car.renew_state(&movement);
        }

        assert_abs_diff_eq!(car.vel[0], 10f64, epsilon = 1e-2)
    }

    #[test]
    fn test_periodic_lane() {
        let length = 100f64;
        let boundary = (
            SimplePlanePair::new(0, [0f64, length]).unwrap(),
            BCspecies::<f64>::Periodic,
        );
        let mut timeiter = ConstStep::<f64>::new(1e-1).unwrap();
        timeiter.set_tmax(40f64).unwrap();
        let mut car = Car::new(0, 1f64, 10f64, 1f64, 0f64);

        for (_t, dt) in timeiter.into_diff() {
            let force = car.drift_force();
            let mut movement = car.euler(&force, dt);
            boundary.check_bc(&car, &mut movement);
            car.renew_state(&movement);
            assert!(car.pos[0] < length);
        }
    }

    #[test]
    fn test_car_sizes() {
        let mut timeiter = ConstStep::<f64>::new(1e-1).unwrap();
        timeiter.set_tmax(100f64).unwrap();
        let mut cars = [
            Car::new(0, 1f64, 10f64, 1f64, 0f64),
            Car::new(0, 1f64, 50f64, 1f64, 0f64),
        ];

        cars[0].pos[0] = 50f64;

        for (_t, dt) in timeiter.into_diff() {
            let mut force1 = cars[0].drift_force();
            let mut force2 = cars[1].drift_force();

            force1 += cars[0].force_from(&cars[1]);
            force2 += cars[1].force_from(&cars[0]);

            let movement1 = cars[0].euler(&force1, dt);
            cars[0].renew_state(&movement1);

            let movement2 = cars[1].euler(&force2, dt);
            cars[1].renew_state(&movement2);

            assert!(cars[0].pos[0] - cars[1].pos[0] > 2.0f64);
            assert!(cars[0].vel[0] < 10f64);
        }

        assert!(cars[1].vel[0] < 10f64);
    }
}