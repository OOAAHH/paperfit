use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write, Read};
use std::path::Path;
use std::collections::HashMap;
use std::sync::Mutex;
use rayon::prelude::*;

/// 计算文件的 MD5 值
fn calculate_md5(file_path: &Path) -> Option<String> {
    let mut file = File::open(file_path).ok()?;
    let mut hasher = md5::Context::new();
    let mut buffer = [0; 4096];

    loop {
        let bytes_read = file.read(&mut buffer).ok()?;
        if bytes_read == 0 {
            break;
        }
        hasher.consume(&buffer[..bytes_read]);
    }

    Some(format!("{:x}", hasher.compute()))
}

/// 读取已有的 MD5 文档并存储为 HashMap
fn load_existing_md5(file_path: &Path) -> std::io::Result<HashMap<String, String>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut md5_map = HashMap::new();

    for line in reader.lines() {
        if let Ok(line) = line {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 2 {
                md5_map.insert(parts[1].to_string(), parts[0].to_string());
            }
        }
    }

    Ok(md5_map)
}

/// 并行计算目录下所有文件的 MD5 并与已有 MD5 值比对
fn calculate_md5_in_directory(
    directory_path: &Path,
    existing_md5_map: &HashMap<String, String>,
    output_mismatch_file: &Path,
    output_unused_file: &Path,
) -> std::io::Result<()> {
    let md5_map = Mutex::new(HashMap::new());

    let paths: Vec<_> = fs::read_dir(directory_path)?
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .collect();

    paths.par_iter().for_each(|file_path| {
        if let Some(md5_value) = calculate_md5(file_path) {
            let relative_path = file_path.strip_prefix(directory_path).unwrap();
            let relative_path_str = relative_path.to_string_lossy().to_string();

            let mut md5_map = md5_map.lock().unwrap();
            md5_map.insert(relative_path_str, md5_value);
        }
    });

    let md5_map = md5_map.into_inner().unwrap();

    let mut mismatch_file = File::create(output_mismatch_file)?;
    for (file_path, md5_value) in &md5_map {
        if let Some(existing_md5_value) = existing_md5_map.get(file_path) {
            if md5_value != existing_md5_value {
                writeln!(mismatch_file, "{}\t{}", md5_value, file_path)?;
            }
        }
    }

    let mut unused_file = File::create(output_unused_file)?;
    for (file_path, existing_md5_value) in existing_md5_map {
        if !md5_map.contains_key(file_path) {
            writeln!(unused_file, "{}\t{}", existing_md5_value, file_path)?;
        }
    }

    Ok(())
}

fn main() -> std::io::Result<()> {
    let directory_path = Path::new("/path/to/your/files");
    let existing_md5_file = Path::new("/path/to/your/existing_md5.txt");
    let output_mismatch_file = Path::new("/path/to/your/mismatch_output.txt");
    let output_unused_file = Path::new("/path/to/your/unused_output.txt");

    let existing_md5_map = load_existing_md5(existing_md5_file)?;

    calculate_md5_in_directory(
        directory_path,
        &existing_md5_map,
        output_mismatch_file,
        output_unused_file,
    )?;

    Ok(())
}
