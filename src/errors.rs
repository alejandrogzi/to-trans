use thiserror;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(thiserror::Error, Debug, PartialEq)]
pub enum Error {
    #[error("Error accessing transcript attribute in: {0}")]
    AttributeError(String),

    #[error("Error parsing attribute: {0}")]
    ParseAttributeError(String),
}

impl From<std::io::Error> for Error {
    fn from(error: std::io::Error) -> Self {
        Error::AttributeError(error.to_string())
    }
}
