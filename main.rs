use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::sync::{Arc, Mutex};
use std::thread;

// 配置线程数量
const NUM_THREADS: usize = 40;

// 定义数据结构
#[derive(Debug, Clone)]
struct PamData {
    sequence: String,
    count: u32,
}

fn main() -> std::io::Result<()> {
    // 获取命令行参数
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <upstream_file> <pam_file>", args[0]);
        std::process::exit(1);
    }
    let upstream_file = &args[1];
    let pam_file = &args[2];
    let output_file = "Fig6E_type1_Results_with_01MatchCount.txt";

    // 读取 upstream 数据
    let upstream_data = read_upstream_data(upstream_file)?;

    // 读取 PAM 数据
    let pam_data = read_pam_data(pam_file)?;

    // 使用 Arc 和 Mutex 来实现线程安全的数据共享
    let pam_data = Arc::new(pam_data);
    let results = Arc::new(Mutex::new(Vec::new()));

    // 创建线程池
    let mut handles = Vec::new();
    let chunk_size = upstream_data.len() / NUM_THREADS + 1;

    for chunk in upstream_data.chunks(chunk_size) {
        let chunk = chunk.to_vec();
        let pam_data = Arc::clone(&pam_data);
        let results = Arc::clone(&results);

        let handle = thread::spawn(move || {
            let mut local_results = Vec::new();
            for data in chunk {
                let (original_line, query_seq) = data;
                let mut perfect_match_count = 0;
                let mut mismatch_1_count = 0;

                // 正向匹配
                for pam in &*pam_data {
                    let mismatches = hamming_distance(&query_seq, &pam.sequence);
                    if mismatches == 0 {
                        perfect_match_count += pam.count;
                    } else if mismatches == 1 {
                        mismatch_1_count += pam.count;
                    }
                }

                // 反向互补匹配
                let rev_comp_seq = reverse_complement(&query_seq);
                for pam in &*pam_data {
                    let mismatches = hamming_distance(&rev_comp_seq, &pam.sequence);
                    if mismatches == 0 {
                        perfect_match_count += pam.count;
                    } else if mismatches == 1 {
                        mismatch_1_count += pam.count;
                    }
                }

                local_results.push(format!(
                    "{}\t{}\t{}",
                    original_line, perfect_match_count, mismatch_1_count
                ));
            }

            // 将结果写入共享的结果集
            let mut results = results.lock().unwrap();
            results.extend(local_results);
        });

        handles.push(handle);
    }

    // 等待所有线程完成
    for handle in handles {
        handle.join().unwrap();
    }

    // 写入结果文件
    let results = results.lock().unwrap();
    let mut file = File::create(output_file)?;
    for line in &*results {
        writeln!(file, "{}", line)?;
    }

    println!("Done! Results saved to {}", output_file);
    Ok(())
}

// 读取 upstream 数据
fn read_upstream_data(filename: &str) -> std::io::Result<Vec<(String, String)>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut data = Vec::new();
    let mut upstream_column_index: Option<usize> = None;

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        // 处理第一行（头部行），检查列名
        if line_num == 0 {
            upstream_column_index = fields.iter().position(|&field| field == "Upstream20bp");
            if upstream_column_index.is_none() {
                eprintln!("Error: 'Upstream20bp' column not found in file!");
                std::process::exit(1);
            }
            continue;
        }

        // 如果找到 Upstream20bp 列名所在的索引，则读取对应列数据
        if let Some(index) = upstream_column_index {
            if fields.len() > index {
                data.push((line.clone(), fields[index].to_string()));
            }
        }
    }
    Ok(data)
}

// 读取 PAM 数据
fn read_pam_data(filename: &str) -> std::io::Result<Vec<PamData>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut data = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() > 2 && fields[0].len() == 20 {
            let sequence = fields[0].to_string();
            let count = fields[2].parse::<u32>().unwrap_or(0);
            data.push(PamData { sequence, count });
        }
    }
    Ok(data)
}

// 计算 Hamming 距离
fn hamming_distance(seq1: &str, seq2: &str) -> usize {
    seq1.chars()
        .zip(seq2.chars())
        .filter(|(a, b)| a != b)
        .count()
}

// 生成反向互补序列的函数
fn reverse_complement(sequence: &str) -> String {
    let mut rev_comp = String::new();
    for base in sequence.chars().rev() {
        rev_comp.push(match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => base, // 处理其他情况（如N或其他字符）
        });
    }
    rev_comp
}
