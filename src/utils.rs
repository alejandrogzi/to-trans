use std::fs::File;
use std::io::{self, Read};
use std::path::PathBuf;

use crate::errors::*;
use crate::gene::GeneModel;

use rayon::prelude::*;

pub fn parallel_parse<'a>(s: &'a str, opt: &str) -> Result<Vec<GeneModel>> {
    let records = s
        .par_lines()
        .filter(|line| !line.starts_with("#"))
        .map(|line| GeneModel::parse(line, opt))
        .collect::<Vec<GeneModel>>();

    return Ok(records);
}

pub fn reader(file: &PathBuf) -> io::Result<String> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}
