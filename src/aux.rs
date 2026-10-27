use std::time::Instant;
use std::fs::{File, OpenOptions};

pub const EMPTY_LEG: i32 = -1;
pub const FLIPPED: i32 = -2;
pub const UNFLIPPED: i32 = -1;

pub const NULL_OP: i32 = -1;
pub const ID_SITE: i32 = -5;
pub const ID_BOND: i32 = -3;    
pub const XX_BOND: i32 = -4;  
pub trait RangeCheck {
    fn within(&self, start: Self, end:Self) -> bool;
}

impl RangeCheck for usize {
    fn within(&self, start: Self, end:Self) -> bool {
        *self >= start && *self < end 
    }
}

#[macro_export]
macro_rules! to_back {
    ($leg: expr) => { $leg ^ 0b10 };
}

#[macro_export]
macro_rules! to_neighbor {
    ($leg: expr) => { $leg ^ 0b01 };
}

#[macro_export]
macro_rules! to_opposite {
    ($leg: expr) => { ($leg ^ 0b01) ^ 0b10 };
}

pub fn report_time(start_time: Instant) {
    let elapsed = start_time.elapsed().as_secs() as u32;
    let hours = elapsed / 3600;
    let minutes = (elapsed % 3600) / 60;
    let seconds = elapsed % 60;
    println!("■ Runtime: {}h-{}min-{}s", hours, minutes, seconds);
}

pub fn create_new_file(file_path: String) -> File {
    OpenOptions::new().write(true).append(true).create(true).open(file_path).unwrap()
}

/* ==================================================================
    These two macros are to drop the safe checks in vectors
================================================================== */
#[macro_export]
macro_rules! unsafe_get {
    ($vec:expr, $idx:expr) => {
        $vec.get_unchecked($idx).clone()
    };

    ($vec:expr, $row:expr, $col:expr) => {
        $vec.get_unchecked($row).get_unchecked($col).clone()
    };
}

#[macro_export]
macro_rules! unsafe_modify {
    ($vec:expr, $idx:expr) => {
        *$vec.get_unchecked_mut($idx)
    };

    ($vec:expr, $row:expr, $col:expr) => {
        *$vec.get_unchecked_mut($row).get_unchecked_mut($col)
    };
}