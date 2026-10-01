# 批量修改 GRASP2018 / grasp-atomic-suite 径向网格

`patch_grasp_grid.py` 可以操作完整原版源码树，或修改本仓库的集中参数文件。
脚本不自动运行编译或提交代码。
脚本只使用 Python 标准库；本工作区仍使用 `graspkit-tools/.venv`。

## 修改文件开头的配置（推荐）

无需通过命令传入数值。在脚本开头的“用户配置”中填写待修改的源码根目录、
核模型参数和运行模式。例如：

```python
GRASP_SOURCE = "/path/to/GRASP2018-grid2990"  # 包含 src/ 的目录
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"  # preview / apply / check
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,
    "n": 2990,
    "h": 0.05,
    "rnt_scale": 2e-6,
    "hp": 0.0,
    "accy": None,
    "point_n": None,
    "point_h": None,
    "point_rnt_scale": None,
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

`None` 表示不修改对应数值，不表示写入 0。以上数值是用法示例，并非适合任何
体系的推荐网格。`GRASP_SOURCE` 默认为 `None`，必须先填入；相对路径以脚本
所在目录为基准，绝对路径和 `~/` 也可以使用。

在服务器上，不带任何参数执行（Bash/fish 均相同）：

```sh
python3 /path/to/grasp-atomic-suite/scripts/patch_grasp_grid.py
```

先设置 `RUN_MODE="preview"` 查看修改；再改成 `"apply"`，执行同一条命令写入；
最后改成 `"check"`，再次执行检查。检查一致时退出 0，仍有待修改文件时退出 1，
配置错误或源码布局不匹配时退出 2。脚本不启动编译。

配置模式下，`BACKUP_DIRECTORY=None` 会将原始文件和校验和保存到
`GRASP_SOURCE/grid-backups/<唯一时间戳>/`，屏幕会打印备份位置。
也可指定一个尚不存在的备份目录。预览和检查不会创建备份；重复应用相同配置
不会创建新备份。恢复时把 `RESTORE_DIRECTORY` 填成打印出的备份目录，
先用 `RUN_MODE="preview"` 查看，再改成 `"apply"` 执行。
恢复模式自动忽略数值配置和 `BACKUP_DIRECTORY`；恢复完成后把
`RESTORE_DIRECTORY` 改回 `None`。恢复后也需要重新编译。

## 多份 GRASP 与 README 的 CMake 编译流程

为每组参数准备一个独立的干净源码副本，例如 `grasp-grid1990/`、
`grasp-grid2990/`。不要复制旧 `build/`、`build-debug/` 缓存或源码内残留的
`.mod/.o` 文件；CMake 缓存记录原目录的绝对路径。原有计算数据单独保存。

可以在每份源码根目录保存一个本脚本副本，设置 `GRASP_SOURCE="."`，
这样每份目录自带对应参数。相对路径与当前 shell 工作目录无关。
原版布局只需脚本本身，不需要整个 suite 或第三方 Python 包。

修改并检查完毕后，按原版 GRASP README 的 CMake 流程手动构建。
先按服务器环境加载 Fortran/MPI/BLAS/LAPACK，再执行：

```sh
cd /path/to/GRASP2018-grid2990
./configure.sh
cd build
make -j4 install
```

这里的 `make` 在 CMake 生成的 `build/` 目录中执行。
`configure.sh` 要求 `build/` 尚不存在；以后同一份源码再次改参数，可使用：

```sh
cd /path/to/GRASP2018-grid2990/build
make clean
make -j4 install
```

全量构建并安装所有实际使用的程序，包括 MPI；配置输出应确认找到 MPI。
每份目录的程序安装到自己的 `bin/`，计算时用对应程序的绝对路径，
例如 `/path/to/GRASP2018-grid2990/bin/rmcdhf_mpi`，避免混用不同配置。

## 原有命令行用法

命令行接口仍保留。**只要提供命令行参数，就完全使用命令行请求，不合入文件
开头的配置**；例如 `--apply` 仍需同时传入路径/参数。以下是兼容用法：

```bash
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py \
  --grasp ../grasp --nnnp 1990 --n 1179 --h 0.025 \
 --rnt-scale 2e-6 --hp 0 --diff
```

修改本仓库时使用显式布局 `--layout atomic-suite`，默认源码路径就是本仓库：

```bash
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py \
  --layout atomic-suite --grasp . --nnnp 2990 --n 1990 --h 0.025 --diff
```

这个布局只修改 `src/lib/libmod/suite_parameters_M.f90`，要求共享初始化模块、
四个 `rmcdhf90*` 入口、RHFS 和 RIS 入口都存在且调用共享初始化；
g_J 复用 RHFS 的网格入口。不会因某个应用缺失而悄悄跳过它。
默认 `--layout grasp2018` 仍要求完整原版源码树，两种布局都会扫描整个 `src/`
拒绝未识别的网格初始化或未同步的容量定义。
本仓库所有有限核程序默认 `FINITE_N=NNNP`，当前为 2990；`NNN1=NNNP+10`
自动派生。修改容量时实际点数自动跟随；显式 `--n` 则将所有程序统一为所给点数。
脚本会拒绝重复容量声明，并依据公共模块的使用情况区分径向参数和局部同名变量。

这是用法示例，不是已验证适合任何元素的推荐网格。
`--nnnp` 设置编译容量，`--n` 设置有限核模型实际点数。
`--rnt-scale 2e-6` 表示 `RNT=2e-6/Z`，不是绝对半径。
未指定的数值保留；有限核的 `N=NNNP` 默认会跟随新的容量。
`--h/--n/--rnt-scale` 只改变有限核分支。
点核分支需要显式使用 `--point-h/--point-n/--point-rnt-scale`；
`--hp` 位于公共赋值处，因此同时影响两种核模型。
`--accy` 可指定独立的数值阈值，否则保留原有公式或固定值（原版为 `H**6`）。
本仓库布局还允许 `--accy 0` 恢复自适应 `H**6`；完整原版布局仍要求正数。

在相同命令后追加 `--apply` 写入；可用 `--backup-dir /path/to/new-directory`
指定一个尚不存在的备份目录，否则备份保存到系统临时目录。
需要长期保存恢复能力时，应明确指定持久化的备份目录。
每次应用保存所有被修改文件的原始字节和校验和。

```bash
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py \
  --grasp ../grasp --nnnp 1990 --n 1179 --h 0.025 \
  --rnt-scale 2e-6 --hp 0 --apply --backup-dir /path/to/new-backup
```

把 `--apply` 换成 `--check` 可以验证补丁已完整应用：一致时退出 0，
仍有待修改文件时退出 1，布局无法识别或参数不合法时退出 2。
重复应用相同配置不会重复插入代码或重新创建备份。

预览恢复，然后追加 `--apply` 执行恢复：

```bash
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py \
  --restore /path/to/new-backup --diff
```

如果源码在应用后又被修改，恢复会拒绝覆盖后续改动。
完整修改计划在写入前生成；所选布局结构不匹配、缺少文件、未处理的径向
容量定义或未知网格初始化会使整个计划失败。
写入使用逐文件原子替换，正常异常会回滚已完成的文件；掉电或强制终止
时可利用备份人工恢复。脚本不支持并发编辑同一个源码树。

原版布局覆盖 17 个径向应用入口、`rwfnrotate/rwfnrelabel` 和
`rwfnestimate` 中重复的容量声明。它还恢复 `rwfnrelabel` 被注释的网格
初始化，并在 RMCDHF 的运行时 ACCY 提示前重新计算默认 ACCY，保留用户
随后显式输入 ACCY 的能力。指定 `--accy` 时也更新 RCI 恢复路径的数值阈值。
不指定 `--rnt-scale` 时，`rwfnrotate` 原有未除以 Z 的 RNT 默认值仍保留。

独立的非相对论 HF、旧 MCHF 绘图/转换工具 `wfnplot/rwfnmchfmcdf`
不在作用范围内。旧 RCI `.res` 恢复时仍读取文件内的网格；要测试新网格应
在新的计算目录启动新计算。交互输入也能覆盖新源码默认值。

本仓库修改时将 `SOURCE_LAYOUT` 设为 `"atomic-suite"`，`GRASP_SOURCE` 指向
本仓库根目录；编译时需开启所有实际使用的串行/MPI 目标。

编译器和 MPI 初始化沿用 GRASP 的构建说明。
确认 PATH 或绝对路径指向新安装的可执行文件，然后检查各阶段 `.sum`
中实际使用的 `RNT/H/HP/N/ACCY`。补丁通过检查仅证明源码配置一致，
不证明实际二进制已更新，也不证明物理量已达到网格收敛。

离线验证脚本：

```bash
../graspkit-tools/.venv/bin/python test/test_patch_grasp_grid.py
```
