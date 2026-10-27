# 强到弱对称性自发破缺

作者：丁一茗

创建日期：2026年6月9日

最后更新：2026年6月9日

## 项目简介

基于量子蒙特卡洛（SSE）模拟**横场伊辛模型（TFIM）**在退相干下的**强-弱对称性自发破缺（SWSSB）**，并计算 Rényi-2 关联函数和对应的 Binder ratio。

论文参考：https://arxiv.org/abs/2603.24342

## 模型/问题

横场伊辛模型：
$$
H = -J \sum_{\langle i,j \rangle} \sigma_i^z \sigma_j^z - h \sum_i \sigma_i^x
$$
在局域 Pauli-Z 噪声（退相干强度 `p`）下，研究：
- **强对称性**：在密度矩阵层面破缺
- **弱对称性**：在期望值层面保持

通过 Rényi-2 关联函数 `<Z_i Z_j>_2` 区分上述两种对称性破缺。

## 算法与方法

基于副本技巧的 SSE QMC 算法，系统总副本数 `Nr = 4*Nb + 2`：

1. **rho0 制备**（2 副本）：标准 SSE QMC 制备热平衡密度矩阵
2. **rho1 构建**（4*Nb 副本）：每个键对应的"小副本"（6 个虚时切片），通过 W 算符编码退相干
3. **核心更新**：
   - **对角更新**：在 rho0 副本中插入/删除对角算符
   - **合并/分裂更新**：按退相干率 `p` 控制副本间的隧道连接
   - **团簇更新**：Swendsen-Wang 型团簇翻转，跨越所有副本
4. **测量**：Rényi-2 ZZ 关联函数、磁化强度、Binder 比

## 依赖环境

```text
Rust（edition 2024，1.85+）
prng_mt 0.1.0（Mersenne Twister PRNG）
Python 3 + NumPy（任务调度和数据收集）
```

安装 Rust：https://rustup.rs

## 文件结构

```text
SWSSB-TFIM/
├── README.md
├── Cargo.toml             # Rust 包配置
├── Cargo.lock             # 依赖锁定
├── .gitignore
├── clear.sh               # 清理所有生成文件
├── run.sh                 # 一键运行（调用 py/run.py）
├── src/                   # Rust 源代码
│   ├── main.rs            # 主程序入口
│   ├── aux.rs             # 辅助常量与宏
│   ├── setup.rs           # SETUP 结构体定义
│   └── setup/
│       ├── init.rs        # 构造函数和初始化
│       ├── measure.rs     # Rényi-2 测量函数
│       ├── random.rs      # 随机数生成
│       ├── qmc.rs         # 更新协调器
│       └── qmc/
│           ├── cluster_update.rs   # 团簇更新
│           ├── diagonal_update.rs  # 对角更新
│           ├── init_reps.rs        # 副本初始化
│           ├── make_rho0.rs        # 单副本热化
│           ├── make_vertices.rs    # 顶点列表构建
│           └── merge_n_split.rs    # 合并/分裂动力学
    └── py/                    # Python 工具
        ├── run.py             # 本地并行运行
        └── collect.py         # 数据收集与合并
```

## 运行方法

### 本地运行
```bash
python3 py/run.py    # 编译并启动并行任务
python3 py/collect.py  # 合并各核心数据
```

### 手动运行
```bash
cargo build --release
./target/release/swssb_1d_unsafe Lx Ly beta J h p N_thm N_rethm N_stat N_bins seed path
```

## 输入参数说明

| 参数 | 描述 |
|------|------|
| `Lx, Ly` | 晶格尺寸（Ly=1 为一维链） |
| `beta` | 逆温度 |
| `J` | 伊辛耦合强度 |
| `h` | 横场强度 |
| `p` | 退相干强度 |
| `N_thm` | 热化步数 |
| `N_rethm` | 副本热化步数 |
| `N_stat` | 每箱测量步数 |
| `N_bins` | 箱数 |
| `seed` | 随机种子 |
| `path` | 输出路径 |

## 输出结果说明

```text
data/<params>/
├── params/               参数文件
└── p<value>/
    ├── core{k}/          各核心原始数据
    │   ├── renyi2_zz.txt     Rényi-2 ZZ 关联函数
    │   ├── renyi2_mag.txt    Rényi-2 磁化强度
    │   ├── renyi2_mag2.txt   Rényi-2 磁化强度平方
    │   ├── renyi2_mag4.txt   Rényi-2 磁化强度四次方
    │   └── renyi2_br.txt     Rényi-2 Binder 比
    └── renyi2_*.txt          收集合并后的数据
logs/                      各核心运行日志
```

## 参考文献

Y.-M. Ding, Y. Guo, Z. Bi, Z. Yan, "Strong-to-Weak Spontaneous Symmetry Breaking in a (2+1)D Transverse-Field Ising Model under Decoherence," arXiv:2603.24342 (2026). https://arxiv.org/abs/2603.24342

## 备注

- Rust edition 2024 需要较新的 rustc 版本。
- 代码中使用了 `unsafe` Rust 以跳过边界检查提升性能。
- `py/run.py` 中可通过修改 `nohup` 标志控制是否后台运行。
