use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write, Read};
use std::path::{Path, PathBuf};
use std::collections::{HashMap, HashSet};
use std::sync::Mutex;
use rayon::prelude::*;
use std::time::Instant; // 用于计时

fn calculate_md5(file_path: &Path) -> Option<(String, u128)> {
    let start_time = Instant::now(); // 开始计时
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

    let md5_value = format!("{:x}", hasher.compute());
    let duration = start_time.elapsed().as_millis(); // 计算耗时（毫秒）

    Some((md5_value, duration))
}

/// 读取已有的 MD5 文档并存储为 HashMap
fn load_existing_md5(file_path: &Path) -> std::io::Result<HashMap<String, (String, String)>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut md5_map = HashMap::new();

    for line in reader.lines() {
        if let Ok(line) = line {
            let parts: Vec<&str> = line.trim().split_whitespace().collect();
            if parts.len() == 2 {
                let md5_value = parts[0].to_string();
                // 去掉前导斜杠，统一路径格式
                let file_path_str = parts[1].trim_start_matches('/').to_string();
                md5_map.insert(file_path_str.clone(), (md5_value, line.clone()));
            }
        }
    }

    Ok(md5_map)
}

/// 获取目录下所有文件的相对路径
fn get_all_files(directory_path: &Path) -> std::io::Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in fs::read_dir(directory_path)? {
        let path = entry?.path();
        if path.is_dir() {
            files.extend(get_all_files(&path)?);
        } else {
            files.push(path);
        }
    }
    Ok(files)
}

/// 主函数，执行 MD5 比对
fn compare_files(
    directory_path: &Path,
    existing_md5_map: &HashMap<String, (String, String)>,
    output_mismatch_file: &Path,
    output_unused_file: &Path,
    output_debug_csv: &Path,
) -> std::io::Result<()> {
    // 获取目录下所有文件的相对路径
    let all_files = get_all_files(directory_path)?;

    let relative_file_paths: Vec<String> = all_files
        .iter()
        .map(|path| {
            // 将文件路径转换为相对于 HRA000087 的相对路径
            // 获取 directory_path 的父目录
            let parent_dir = directory_path.parent().unwrap();
            // 从父目录开始计算相对路径
            let rel_path = path.strip_prefix(parent_dir).unwrap();
            // 转换为字符串并去掉可能的前导斜杠
            rel_path.to_string_lossy().trim_start_matches('/').to_string()
        })
        .collect();

    // 将相对路径放入 HashSet，方便查找
    let file_set: HashSet<String> = relative_file_paths.iter().cloned().collect();

    // 找出在目录中存在但不在 md5sum.txt 中的文件，记录到 unused_output.txt
    {
        let mut unused_file = File::create(output_unused_file)?;
        for file_path in &relative_file_paths {
            if !existing_md5_map.contains_key(file_path) {
                writeln!(unused_file, "{}", file_path)?;
            }
        }
    }

    // 准备调试信息 CSV 文件
    let debug_csv = Mutex::new(csv::Writer::from_path(output_debug_csv)?);

    // 找出在目录和 md5sum.txt 中都存在的文件，进行 MD5 比对
    {
        let mismatch_file = Mutex::new(File::create(output_mismatch_file)?);

        // 创建一个迭代器，包含需要比对的文件及其相对路径
        let files_to_check: Vec<_> = all_files
            .iter()
            .zip(relative_file_paths.iter())
            .filter(|(_, rel_path)| existing_md5_map.contains_key(*rel_path))
            .collect();

        // 并行计算 MD5 并比对
        files_to_check.par_iter().for_each(|(file_path, rel_path)| {
            if let Some((md5_value, duration)) = calculate_md5(file_path) {
                if let Some((existing_md5_value, md5_line)) = existing_md5_map.get(*rel_path) {
                    // 写入调试信息 CSV 文件
                    let mut debug_csv = debug_csv.lock().unwrap();
                    let record = DebugInfo {
                        computed_md5: md5_value.clone(),
                        file_path: (*rel_path).clone(),
                        md5_line: md5_line.clone(),
                        duration,
                    };
                    debug_csv.serialize(record).unwrap();

                    if &md5_value != existing_md5_value {
                        let mut mismatch_file = mismatch_file.lock().unwrap();
                        writeln!(mismatch_file, "{}\t{}", md5_value, rel_path).unwrap();
                        println!("文件 {} 的 MD5 值不一致", rel_path);
                    }
                }
            }
        });
    }

    Ok(())
}

#[derive(serde::Serialize)]
struct DebugInfo {
    computed_md5: String,
    file_path: String,
    md5_line: String,
    duration: u128,
}

fn main() -> std::io::Result<()> {
    let directory_path = Path::new("/hdd/HRA000087/HRA000087/");
    let existing_md5_file = directory_path.join("md5sum.txt");
    let output_mismatch_file = directory_path.join("mismatch_output.txt");
    let output_unused_file = directory_path.join("unused_output.txt");
    let output_debug_csv = directory_path.join("debug_info.csv");

    let existing_md5_map = load_existing_md5(&existing_md5_file)?;

    compare_files(
        directory_path,
        &existing_md5_map,
        &output_mismatch_file,
        &output_unused_file,
        &output_debug_csv,
    )?;

    Ok(())
}
