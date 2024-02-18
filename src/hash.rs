use tiny_keccak::Hasher;

pub fn keccak(data: &[u8]) -> [u8; 32] {
    let mut hasher = tiny_keccak::Keccak::v256();
    let mut result = [0u8; 32];
    hasher.update(data);
    hasher.finalize(&mut result);
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_keccak() {
        let result = keccak(b"");
        let result_hex = hex::encode(result);
        assert_eq!(
            result_hex,
            "c5d2460186f7233c927e7db2dcc703c0e500b653ca82273b7bfad8045d85a470"
        );
    }
}
