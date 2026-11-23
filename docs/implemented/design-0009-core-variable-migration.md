# libcore 变量与常量命名迁移

下表记录这批迁移当时的名称。后续改动已将重复的共享状态别名移至定义端，并更新了库过程名。[后续实施与回归](../reviews/evidence-0015-fortran-alias-audit.md#清理后的命名与验证) · [现行过程对照](design-0007-shared-library-procedure-names.md)

日期 2026-11-20。源码基准为 `2d4c721` 加同一工作树中的 CSF 查询模块合并。
逐文件盘点 CMake 清单中的 71 份 Fortran 源码，覆盖 1,181 个作用域内的声明实体。改动覆盖 68 份源码，迁移了 885 个原名（含 247 个形参、14 个 `PARAMETER` 常量），另有 100 处共享符号引用使用了更明确的局部名称。

模块和源码的原名见[公共基础库命名迁移](design-0008-core-library-naming.md)。

## 阅读方式

- 每个文件下按模块、过程或类型作用域列出。表中「形参」可用于 Fortran 关键字实参。采用关键字调用时，把旧关键字替换为表中现名。位置调用的参数顺序和类型沿用原接口。
- 命名区分轨道与 CSF 索引、2J 与 2J+1、张量阶数和记录位置，径向例程分别使用 `radial_index`、`radius`、`grid_point_count` 等角色名称。物理量的单位已由代码确定时保留单位后缀，例如 `fermi_diffuseness_fm`。
- 逐作用域列出的保留项包括已有明确语义的标识符、公开类型的成员、沿用过程名的隐式返回值，以及只见于声明的历史变量。后者标记为「仅声明」。`npol.inp` 的 `k1`、`b1`、`k2`、`b2` 等键与 NAMELIST 变量名绑定，沿用文件字段名。
- `USE ONLY` 别名的「定义」列标出共享状态的实际拥有模块。代码依然读写该模块的同一对象。Fortran 过程名、C 链接名、I/O 字符串和文件字段沿用已有契约。

## 验证

- 使用 GNU Fortran 15 对比 71 份源码新旧版本解析后的可执行语句。还原标识符后逐文件一致。另核对了声明类型、形参顺序、别名绑定和字符串字面值。
- 启用 MPI 和 OpenMP 的 Release 构建通过 110 项测试。最后一轮名称校正后重新构建，其中 12 项相关 Release 回归重跑通过。关闭 MPI 和 OpenMP 的 Check 构建通过 70 项测试，Rust 读库测试在 Release 配置中执行。

## [angular_coefficient_buffer.f90](../../src/lib/core/angular_coefficient_buffer.f90)

### `angular_coefficient_buffer → alcbuf`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `mode` | `buffer_operation` | 形参 |
| `newsiz` | `new_buffer_capacity` | 变量 |

保留 `ierr`（仅声明）。

## [nonnegative_atan2.f90](../../src/lib/core/nonnegative_atan2.f90)

### `nonnegative_atan2 → arctan`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `arg1` | `y_coordinate` | 形参 |
| `arg2` | `x_coordinate` | 形参 |
| `first` | `check_intrinsic_range` | 变量 |
| `intrin` | `intrinsic_returns_nonnegative_angle` | 变量 |

保留 `arctan`。

## [calendar_clock.f90](../../src/lib/core/calendar_clock.f90)

### `calendar_clock → calen`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `jtime` | `time_text` | 形参 |
| `jdate` | `date_text` | 形参 |

## [complex_gamma.f90](../../src/lib/core/complex_gamma.f90)

### `complex_gamma → cgamma`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `argi` | `argument_imaginary` | 形参 |
| `argr` | `argument_real` | 形参 |
| `resr` | `gamma_real` | 形参 |
| `resi` | `gamma_imaginary` | 形参 |
| `i` | `stirling_term_index` | 变量 |
| `fn` | `stirling_coefficients` | 变量 |
| `fd` | `bernoulli_denominators` | 变量 |
| `hlntpi` | `half_log_two_pi` | 变量 |
| `twoi` | `twice_term_index` | 变量 |
| `diff` | `distance_to_integer` | 变量 |
| `argum` | `argument_magnitude` | 变量 |
| `clngi` | `log_gamma_imaginary` | 变量 |
| `facneg` | `log_reflection_factor` | 变量 |
| `argur` | `shifted_argument_real` | 变量 |
| `ovlfac` | `recurrence_product` | 变量 |
| `clngr` | `log_gamma_real` | 变量 |
| `fac` | `scale_factor` | 变量 |
| `obasq` | `inverse_argument_squared` | 变量 |
| `argui` | `shifted_argument_imaginary` | 变量 |
| `argui2` | `imaginary_part_squared` | 变量 |
| `ovlfr` | `recurrence_product_real` | 变量 |
| `ovlfi` | `recurrence_product_imaginary` | 变量 |
| `termr` | `term_real` | 变量 |
| `termi` | `term_imaginary` | 变量 |
| `argur2` | `real_part_squared` | 变量 |
| `obasqr` | `inverse_square_real` | 变量 |
| `obasqi` | `inverse_square_imaginary` | 变量 |
| `zfacr` | `argument_power_real` | 变量 |
| `zfaci` | `argument_power_imaginary` | 变量 |
| `first` | `initialize_stirling_coefficients` | 变量 |
| `negarg` | `use_reflection_formula` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `cgamma` | `def_c::expmax` | `expmax` | `maximum_exp_argument` |
| `cgamma` | `def_c::expmin` | `expmin` | `minimum_exp_argument` |
| `cgamma` | `def_c::precis` | `precis` | `machine_precision` |

## [relativistic_three_j.f90](../../src/lib/core/relativistic_three_j.f90)

### `relativistic_three_j → clrx`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `kappaa` | `kappa_a` | 形参 |
| `kappab` | `kappa_b` | 形参 |
| `k` | `tensor_rank` | 形参 |
| `ka` | `abs_kappa_a` | 变量 |
| `kb` | `abs_kappa_b` | 变量 |
| `kapkb` | `abs_kappa_sum` | 变量 |
| `kabkp` | `abs_kappa_sum_plus_rank` | 变量 |
| `kamkb` | `abs_kappa_a_minus_b` | 变量 |
| `kbmka` | `abs_kappa_b_minus_a` | 变量 |
| `exptrm` | `log_magnitude` | 变量 |

保留 `clrx`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `clrx` | `facts_c::gam` | `gam` | `log_factorials` |

## [integer_format.f90](../../src/lib/core/integer_format.f90)

### `integer_format → convrt`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `intnum` | `integer_value` | 形参 |
| `lenth` | `text_length` | 形参 |
| `cnum` | `number_text` | 形参 |
| `form` | `integer_format` | 变量 |
| `c1020` | `two_digit_widths` | 变量 |
| `c19` | `one_digit_widths` | 变量 |

## [integer_format_with_context.f90](../../src/lib/core/integer_format_with_context.f90)

### `integer_format_with_context → convrt2`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `intnum` | `integer_value` | 形参 |
| `lenth` | `text_length` | 形参 |
| `cnum` | `number_text` | 形参 |
| `from` | `caller_name` | 形参 |
| `form` | `integer_format` | 变量 |
| `c1020` | `two_digit_widths` | 变量 |
| `c19` | `one_digit_widths` | 变量 |

## [half_integer_format.f90](../../src/lib/core/half_integer_format.f90)

### `half_integer_format → convrt_double`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `intnum` | `twice_angular_momentum` | 形参 |
| `lenth` | `text_length` | 形参 |
| `cnum` | `angular_momentum_text` | 形参 |
| `intnumgg` | `numerator_to_print` | 变量 |
| `form` | `integer_format` | 变量 |
| `c1020` | `two_digit_widths` | 变量 |
| `c19` | `one_digit_widths` | 变量 |

## [closed_shell_mcp.f90](../../src/lib/core/closed_shell_mcp.f90)

### `closed_shell_mcp → cord`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ja` | `bra_csf_index` | 形参 |
| `jb` | `ket_csf_index` | 形参 |
| `ja1` | `bra_shell_position` | 形参 |
| `jb1` | `ket_shell_position` | 形参 |
| `ipca` | `orbital_list_selector` | 形参 |
| `eps` | `coefficient_cutoff` | 常量 |
| `ia1` | `first_orbital` | 变量 |
| `ib1` | `second_orbital` | 变量 |
| `ns` | `swap_orbital` | 变量 |
| `kap1` | `first_kappa` | 变量 |
| `j1` | `first_abs_kappa` | 变量 |
| `nqs1` | `first_occupation` | 变量 |
| `numax` | `maximum_rank` | 变量 |
| `nu` | `tensor_rank` | 变量 |
| `kap2` | `second_kappa` | 变量 |
| `j2` | `second_abs_kappa` | 变量 |
| `nqs2` | `second_occupation` | 变量 |
| `numin` | `minimum_rank` | 变量 |
| `x` | `mcp_coefficient` | 变量 |
| `const` | `occupation_factor` | 变量 |
| `gam` | `three_j_symbol` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `cord` | `m_c::jlist` | `jlist` | `peel_orbital_indices` |
| `cord` | `m_c::klist` | `klist` | `core_orbital_indices` |
| `cord` | `m_c::nq1` | `nq1` | `bra_occupations` |

## [radial_node_count.f90](../../src/lib/core/radial_node_count.f90)

### `radial_node_count → count_radial_nodes`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `fr` | `radial_values` | 形参 |
| `mtpfr` | `last_radial_index` | 形参 |
| `nncff` | `node_count` | 形参 |
| `sgn` | `first_extremum_sign` | 形参 |
| `lcext` | `extremum_grid_indices` | 变量 |
| `next` | `extremum_count` | 变量 |
| `i` | `sample_index` | 变量 |
| `loc` | `extremum_index` | 变量 |
| `nstps` | `remaining_extrema` | 变量 |
| `trace_env` | `trace_setting` | 变量 |
| `trace_count` | `trace_enabled` | 变量 |
| `trace_init` | `trace_initialized` | 变量 |
| `ext` | `segment_max_amplitude` | 变量 |
| `emx` | `peak_amplitude` | 变量 |
| `abfri` | `sample_amplitude` | 变量 |
| `test` | `absolute_sign_sum` | 变量 |
| `threse` | `extremum_cutoff` | 变量 |
| `ablcl` | `extremum_amplitude` | 变量 |

保留 `trace_status`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `count_impl` | `coun_c::thresh` | `thresh` | `relative_extremum_cutoff` |
| `count_impl` | `def_c::accy` | `accy` | `working_tolerance` |

## [reduced_angular_matrix_element.f90](../../src/lib/core/reduced_angular_matrix_element.f90)

### `reduced_angular_matrix_element → cre`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `kap1` | `bra_kappa` | 形参 |
| `kap2` | `ket_kappa` | 形参 |
| `k` | `tensor_rank` | 形参 |
| `k1` | `abs_bra_kappa` | 变量 |
| `dk1k2` | `subshell_degeneracy_product` | 变量 |

保留 `cre`。

## [csf_full_basis_reader.f90](../../src/lib/core/csf_full_basis_reader.f90)

### `csf_full_basis_reader → load_csf_basis`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `name` | `state_name` | 形参 |
| `ncore` | `core_subshell_count` | 形参 |
| `source` | `basis_reader` | 变量 |
| `b` | `block_index` | 变量 |
| `first` | `first_csf` | 变量 |
| `last` | `last_csf` | 变量 |
| `i` | `source_orbital_index` | 变量 |
| `kept` | `retained_orbital_count` | 变量 |
| `offset` | `csf_offset` | 变量 |
| `new_quantum` | `new_subshell_quantum_codes` | 变量 |

保留 `status`、`new_occupations`、`new_couplings`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `csf_basis_input` | `blk_c::nblock` | `nblock` | `symmetry_block_count` |
| `csf_basis_input` | `blk_c::ncfblk` | `ncfblk` | `cumulative_csf_counts` |

## [csf_basis_merging.f90](../../src/lib/core/csf_basis_merging.f90)

### `csf_basis_merging → merge_csf_bases`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `left_core` | `left_core_subshell_count` | 形参 |
| `right_core` | `right_core_subshell_count` | 形参 |
| `merged_n` | `merged_principal_n` | 变量 |
| `merged_k` | `merged_kappa` | 变量 |
| `left` | `left_orbital_index` | 变量 |
| `right` | `right_orbital_index` | 变量 |
| `choice` | `merge_choice` | 变量 |
| `count` | `merged_orbital_count` | 变量 |
| `i` | `merged_orbital_index` | 变量 |
| `j` | `input_orbital_index` | 变量 |
| `rows` | `merged_csf_count` | 变量 |
| `wide_rows` | `wide_merged_csf_count` | 变量 |
| `left_common` | `left_orbital_in_right_tail` | 变量 |
| `right_common` | `right_orbital_in_left_tail` | 变量 |
| `quantum` | `subshell_quantum_codes` | 变量 |

保留 `status`、`message`、`merged_labels`、`allocation_status`、`occupations`、`couplings`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `csf_basis_merge` | `def_c::nelecr` | `nelecr` | `left_electron_count` |
| `csf_basis_merge` | `orb_c::ncfr` | `ncfr` | `left_csf_count` |
| `csf_basis_merge` | `orb_c::nwr` | `nwr` | `left_orbital_count` |
| `csf_basis_merge` | `orb_c::npr` | `npr` | `left_principal_n` |
| `csf_basis_merge` | `orb_c::nakr` | `nakr` | `left_kappa` |
| `csf_basis_merge` | `orb_c::nhr` | `nhr` | `left_orbital_labels` |
| `csf_basis_merge` | `orb_c::iqar` | `iqar` | `left_occupations` |
| `csf_basis_merge` | `stat_c::jqsar` | `jqsar` | `left_subshell_quantum_codes` |
| `csf_basis_merge` | `stat_c::jcupar` | `jcupar` | `left_coupling_codes` |

## [rci_csf_block_reader.f90](../../src/lib/core/rci_csf_block_reader.f90)

### `rci_csf_block_reader`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `source` | `block_reader` | 变量 |

### `rci_csf_block_reader → open_csf_blocks`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `name` | `state_name` | 形参 |
| `ncore` | `core_subshell_count` | 形参 |
| `labels` | `block_labels` | 形参 |

保留 `status`。

### `rci_csf_block_reader → read_csf_zero_order_counts`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `counts` | `zero_order_counts` | 形参 |

保留 `status`。

### `rci_csf_block_reader → load_csf_block`

保留 `block_index`、`status`。

## [csf_descriptor_reader.f90](../../src/lib/core/csf_descriptor_reader.f90)

### `csf_descriptor_reader → descriptor_source`

保留 `handle`、`is_open`、`header_loaded`、`core_count`、`peel_count`、`orbital_count`、`block_count`、`row_width`、`total_csf_count`、`electron_count`、`next_block_index`、`principal_numbers`、`relativistic_kappas`、`orbital_momenta`、`single_electron_two_j`、`orbital_labels`、`block_lengths`、`block_two_j`、`block_parities`。

### `csf_descriptor_reader → failed`

保留 `failed`、`status`。

### `csf_descriptor_reader → require_descriptor`

保留 `ok`、`condition`、`message`、`status`。

### `csf_descriptor_reader → open_descriptor_source`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `result` | `open_status` | 变量 |

保留 `source`、`filename`、`maximum_blocks`、`header_filename`、`status`。

### `csf_descriptor_reader → open_descriptor_source_impl`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `header` | `header_path` | 变量 |
| `stem` | `descriptor_stem` | 变量 |
| `message` | `reader_message` | 变量 |
| `info` | `header_counts` | 变量 |
| `lengths` | `block_lengths` | 变量 |
| `orbitals` | `orbital_table` | 变量 |
| `symmetries` | `symmetry_table` | 变量 |
| `i` | `orbital_index` | 变量 |
| `k` | `kappa` | 变量 |
| `l` | `orbital_l_value` | 变量 |
| `n` | `path_length` | 变量 |
| `letters` | `orbital_letters` | 常量 |

保留 `source`、`filename`、`maximum_blocks`、`header_filename`、`status`、`read_status`、`allocation_status`。

### `csf_descriptor_reader → read_descriptor_zero_order_counts`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `values` | `raw_zero_order_counts` | 变量 |
| `message` | `reader_message` | 变量 |

保留 `source`、`counts`、`status`、`read_status`。

### `csf_descriptor_reader → close_descriptor_source`

保留 `source`。

### `csf_descriptor_reader → descriptor_block_label`

保留 `total_two_j`、`parity`、`label`、`angular_text`、`padded_text`、`parity_symbol`。

### `csf_descriptor_reader → install_descriptor_header`

保留 `source`、`core_count`、`block_labels`、`status`、`block_index`、`allocation_status`、`prepared_lengths`。

### `csf_descriptor_reader → install_descriptor_header → resolve_subshell_state`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `capacity` | `subshell_capacity` | 变量 |

保留 `single_two_j`、`occupation`、`physical_seniority`、`printed_two_j`、`block_index`、`row_index`、`status`、`subshell_two_j`、`population`、`table_row`、`table_entry`、`table_offset`、`match_count`、`matching_j_count`、`candidate_two_j`、`candidate_seniority`、`candidate_label`。

### `csf_descriptor_reader → install_descriptor_header → resolve_subshell_state → require_state`

保留 `ok`、`condition`、`message`。

### `csf_descriptor_reader → install_descriptor_header → decode_descriptor_row`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `row` | `packed_row` | 形参 |
| `capacity` | `subshell_capacity` | 变量 |

保留 `source`、`block_index`、`row_index`、`occupations`、`subshell_quantum_numbers`、`couplings`、`status`、`shell_two_j`、`printed_two_k`、`peel_index`、`orbital_index`、`channel_offset`、`occupation`、`seniority`、`printed_two_j`、`occupied_count`、`occupied_position`、`open_count`、`previous_two_j`、`coupled_two_j`、`last_nonzero_shell`、`parity_exponent`、`total_two_j`、`parity`、`electron_count`。

### `csf_descriptor_reader → install_descriptor_header → decode_descriptor_row → require_row`

保留 `ok`、`condition`、`message`。

### `csf_descriptor_reader → install_descriptor_header → read_descriptor_block`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `result` | `block_status` | 变量 |

保留 `source`、`block_index`、`occupations`、`subshell_quantum_numbers`、`couplings`、`status`。

### `csf_descriptor_reader → install_descriptor_header → read_descriptor_block_impl`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `row_quantum` | `row_subshell_quantum_numbers` | 变量 |
| `message` | `reader_message` | 变量 |
| `batch_count` | `batch_row_count` | 变量 |
| `batch_row` | `batch_row_index` | 变量 |

保留 `source`、`block_index`、`occupations`、`subshell_quantum_numbers`、`couplings`、`row_couplings`、`status`、`batch_values`、`batch_capacity`、`row_offset`、`read_status`、`allocation_status`。

## [csf_quantum_queries.f90](../../src/lib/core/csf_quantum_queries.f90)

### `csf_quantum_queries → iq`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `isubsh` | `subshell_index` | 形参 |
| `icsf` | `csf_index` | 形参 |

保留 `iq`。

### `csf_quantum_queries → jqs`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `iwhich` | `quantum_number_selector` | 形参 |
| `isubsh` | `subshell_index` | 形参 |
| `icsf` | `csf_index` | 形参 |

保留 `jqs`。

### `csf_quantum_queries → jcup`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `loc` | `coupling_index` | 形参 |
| `icsf` | `csf_index` | 形参 |

保留 `jcup`。

### `csf_quantum_queries → ispar`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `icsf` | `csf_index` | 形参 |

保留 `ispar`。

### `csf_quantum_queries → itjpo`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `icsf` | `csf_index` | 形参 |

保留 `itjpo`。

## [csf_descriptor_c_bindings.f90](../../src/lib/core/csf_descriptor_c_bindings.f90)

### `csf_descriptor_c_bindings → rust_descriptor_version`

保留 `status`。

### `csf_descriptor_c_bindings → rust_descriptor_read_radial_header`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `header` | `header_path` | 形参 |
| `orbitals` | `orbital_table` | 形参 |
| `info` | `header_counts` | 形参 |
| `message` | `error_message` | 形参 |
| `capacity` | `message_capacity` | 形参 |

保留 `orbital_capacity`、`status`。

### `csf_descriptor_c_bindings → rust_descriptor_open`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `path` | `descriptor_path` | 形参 |
| `header` | `header_path` | 形参 |
| `message` | `error_message` | 形参 |
| `capacity` | `message_capacity` | 形参 |

保留 `handle`、`status`。

### `csf_descriptor_c_bindings → rust_descriptor_info`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `info` | `header_counts` | 形参 |
| `message` | `error_message` | 形参 |
| `capacity` | `message_capacity` | 形参 |

保留 `handle`、`status`。

### `csf_descriptor_c_bindings → rust_descriptor_tables`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `orbitals` | `orbital_table` | 形参 |
| `symmetries` | `symmetry_table` | 形参 |
| `lengths` | `block_lengths` | 形参 |
| `message` | `error_message` | 形参 |
| `capacity` | `message_capacity` | 形参 |

保留 `handle`、`orbital_capacity`、`block_capacity`、`status`。

### `csf_descriptor_c_bindings → rust_descriptor_zero_order_counts`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `counts` | `zero_order_counts` | 形参 |
| `message` | `error_message` | 形参 |
| `capacity` | `message_capacity` | 形参 |

保留 `handle`、`block_capacity`、`status`。

### `csf_descriptor_c_bindings → rust_descriptor_read`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `values` | `row_values` | 形参 |
| `size` | `value_capacity` | 形参 |
| `rows` | `row_count` | 形参 |
| `capacity` | `message_capacity` | 形参 |
| `message` | `error_message` | 形参 |

保留 `handle`、`status`。

### `csf_descriptor_c_bindings → rust_descriptor_close`

保留 `handle`。

## [dirac_coulomb_radial_orbital.f90](../../src/lib/core/dirac_coulomb_radial_orbital.f90)

### `dirac_coulomb_radial_orbital → dcbsrw`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `n` | `principal_n` | 形参 |
| `mtp` | `last_radial_index` | 形参 |
| `z` | `effective_nuclear_charge` | 形参 |
| `e` | `binding_energy` | 形参 |
| `rg0` | `large_component_origin_coefficient` | 形参 |
| `rg` | `large_component` | 形参 |
| `rf` | `small_component` | 形参 |
| `k` | `abs_kappa`、`series_index` | 变量 |
| `nr` | `radial_quantum_number` | 变量 |
| `nrfac` | `radial_quantum_factorial` | 变量 |
| `i` | `sample_index` | 变量 |
| `iordr1` | `first_series_order` | 变量 |
| `iordr2` | `second_series_order` | 变量 |
| `alfa` | `fine_structure_constant` | 变量 |
| `fn` | `principal_n_real` | 变量 |
| `fkappa` | `kappa_real` | 变量 |
| `fk` | `abs_kappa_real` | 变量 |
| `fnr` | `radial_quantum_number_real` | 变量 |
| `zalfa` | `charge_times_alpha` | 变量 |
| `gamma` | `origin_exponent` | 变量 |
| `twogp1` | `twice_exponent_plus_one` | 变量 |
| `bign` | `effective_principal_n` | 变量 |
| `eps` | `relativistic_energy_ratio` | 变量 |
| `argr` | `gamma_argument_real` | 变量 |
| `argi` | `gamma_argument_imaginary` | 变量 |
| `rgamm1` | `first_gamma_value` | 变量 |
| `dummy` | `unused_gamma_imaginary` | 变量 |
| `rgamm2` | `second_gamma_value` | 变量 |
| `fac` | `scale_factor` | 变量 |
| `fg` | `large_component_scale` | 变量 |
| `ff` | `small_component_scale` | 变量 |
| `facn` | `term_factorial` | 变量 |
| `a` | `numerator_parameter` | 变量 |
| `an1` | `first_numerator_product` | 变量 |
| `an2` | `second_numerator_product` | 变量 |
| `b` | `denominator_parameter` | 变量 |
| `bn` | `denominator_product` | 变量 |
| `fden` | `inverse_series_denominator` | 变量 |
| `bignmk` | `effective_n_minus_kappa` | 变量 |
| `rho` | `scaled_radius` | 变量 |
| `rhon` | `radius_power` | 变量 |
| `f1` | `first_radial_series` | 变量 |
| `f2` | `second_radial_series` | 变量 |
| `ovlfac` | `radial_envelope` | 变量 |
| `cutoff` | `amplitude_cutoff` | 变量 |

保留 `kappa`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `dcbsrw` | `grid_c::radial_point_count` | `ntp` | `radial_point_count` |
| `dcbsrw` | `tatb_c::ta` | `ta` | `first_series_coefficients` |
| `dcbsrw` | `tatb_c::tb` | `tb` | `second_series_coefficients` |

## [real_vector_initialization.f90](../../src/lib/core/real_vector_initialization.f90)

### `real_vector_initialization → dinit`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `n` | `element_count` | 形参 |
| `incx` | `stride` | 形参 |
| `a` | `fill_value` | 形参 |
| `x` | `values` | 形参 |
| `xaddr` | `first_index` | 变量 |

保留 `i`（仅声明）。

## [sparse_dot_and_scatter.f90](../../src/lib/core/sparse_dot_and_scatter.f90)

### `sparse_dot_and_scatter → dmerge`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `n` | `element_count` | 形参 |
| `dconst` | `update_scale` | 形参 |
| `dl` | `dot_product` | 形参 |
| `idy` | `target_indices` | 形参 |
| `db` | `source_vector` | 形参 |
| `dc` | `updated_vector` | 形参 |
| `da` | `matrix_values` | 形参 |
| `i` | `entry_index` | 变量 |
| `dsum` | `dot_product_sum` | 变量 |

## [radial_wavefunction_derivatives.f90](../../src/lib/core/radial_wavefunction_derivatives.f90)

### `radial_wavefunction_derivatives → dpbdt`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `j` | `orbital_index` | 形参 |
| `i` | `radial_index` | 变量 |
| `irow` | `stencil_row` | 变量 |
| `a1` | `offset_6_weight` | 变量 |
| `a2` | `offset_5_weight` | 变量 |
| `a3` | `offset_4_weight` | 变量 |
| `a4` | `offset_3_weight` | 变量 |
| `a5` | `offset_2_weight` | 变量 |
| `a6` | `offset_1_weight` | 变量 |
| `hdpbdt` | `scaled_large_derivative` | 变量 |
| `hdqbdt` | `scaled_small_derivative` | 变量 |

保留 `k`（仅声明）、`loc`（仅声明）、`aik`（仅声明）。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `dpbdt` | `lic13_c::a` | `a` | `derivative_weights` |
| `dpbdt` | `tatb_c::ta` | `ta` | `scaled_large_derivative_grid` |
| `dpbdt` | `tatb_c::tb` | `tb` | `scaled_small_derivative_grid` |

## [radial_ascii_plot.f90](../../src/lib/core/radial_ascii_plot.f90)

### `radial_ascii_plot → draw`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `mf` | `last_radial_index` | 形参 |
| `sp` | `large_component_scale` | 形参 |
| `sq` | `small_component_scale` | 形参 |
| `p` | `large_component` | 形参 |
| `q` | `small_component` | 形参 |
| `i` | `plot_index` | 变量 |
| `ioffst` | `zero_baseline_row` | 变量 |
| `ixlocq` | `small_component_row` | 变量 |
| `ixlocp` | `large_component_row` | 变量 |
| `iyloc` | `radius_column` | 变量 |
| `fift4` | `plot_height_scale` | 变量 |
| `ohtht` | `plot_width_scale` | 变量 |
| `dmx` | `maximum_amplitude` | 变量 |
| `dmn` | `minimum_amplitude` | 变量 |
| `spi` | `scaled_large_value` | 变量 |
| `sqi` | `scaled_small_value` | 变量 |
| `rmx` | `maximum_radius` | 变量 |
| `xscal` | `amplitude_to_row` | 变量 |
| `yscal` | `index_to_column` | 变量 |
| `sqx` | `small_component_to_row` | 变量 |
| `spx` | `large_component_to_row` | 变量 |
| `first` | `initialize_plot_templates` | 变量 |
| `cblank` | `bordered_row` | 变量 |
| `cdash` | `dashed_row` | 变量 |
| `cplot` | `plot_rows` | 变量 |
| `cp` | `large_component_marker` | 变量 |
| `cq` | `small_component_marker` | 变量 |

保留 `sixty`。

## [alternating_exponential_series.f90](../../src/lib/core/alternating_exponential_series.f90)

### `alternating_exponential_series → es`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `f` | `exponent_argument` | 形参 |
| `s2f` | `inverse_square_sum` | 形参 |
| `s3f` | `inverse_cube_sum` | 形参 |
| `n` | `term_index` | 变量 |
| `fase` | `alternating_sign` | 变量 |
| `en` | `term_index_real` | 变量 |
| `obn` | `inverse_term_index` | 变量 |
| `enf` | `exponential_weight` | 变量 |
| `term2` | `inverse_square_term` | 变量 |
| `term3` | `inverse_cube_term` | 变量 |
| `s2last` | `previous_square_sum` | 变量 |

## [log_factorial_table.f90](../../src/lib/core/log_factorial_table.f90)

### `log_factorial_table → factt`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `i` | `factorial_table_index` | 变量 |
| `x` | `factorial_multiplier` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `factt` | `facts_c::gam` | `gam` | `log_factorials` |
| `factt` | `facts_c::mfact` | `mfact` | `log_factorial_table_size` |

## [subshell_label_selection.f90](../../src/lib/core/subshell_label_selection.f90)

### `subshell_label_selection → getrsl`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `nsubs` | `matched_subshell_count` | 形参 |
| `indx` | `matched_orbital_indices` | 形参 |
| `istart` | `token_start` | 变量 |
| `i` | `character_index` | 变量 |
| `iend` | `token_end` | 变量 |
| `nform` | `principal_digit_count` | 变量 |
| `j` | `orbital_index` | 变量 |
| `lenth` | `width_text_length` | 变量 |
| `ios` | `read_status` | 变量 |
| `npq` | `requested_principal_n` | 变量 |
| `k` | `match_index` | 变量 |
| `found` | `subshell_matched` | 变量 |
| `nlany` | `any_orbital_l` | 变量 |
| `nlok` | `orbital_l_matches` | 变量 |
| `npany` | `any_principal_n` | 变量 |
| `npok` | `principal_n_matches` | 变量 |
| `nsany` | `any_spin_component` | 变量 |
| `nsok` | `spin_component_matches` | 变量 |
| `nlq` | `requested_orbital_letter` | 变量 |
| `nsq` | `requested_spin_suffix` | 变量 |
| `reci` | `current_character` | 变量 |
| `cnum` | `width_text` | 变量 |
| `form` | `integer_format` | 变量 |
| `record` | `input_line` | 变量 |

## [yes_no_input.f90](../../src/lib/core/yes_no_input.f90)

### `yes_no_input → getyn`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `rspns` | `response_character` | 变量 |

保留 `getyn`。

## [one_electron_csf_filter.f90](../../src/lib/core/one_electron_csf_filter.f90)

### `one_electron_csf_filter → ichkq1`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ja` | `bra_csf_index` | 形参 |
| `jb` | `ket_csf_index` | 形参 |
| `k` | `differing_subshell_count` | 变量 |
| `i` | `orbital_index` | 变量 |
| `i_iqa` | `bra_occupation` | 变量 |
| `i_iqb` | `ket_occupation` | 变量 |

保留 `ichkq1`。

## [two_electron_csf_filter.f90](../../src/lib/core/two_electron_csf_filter.f90)

### `two_electron_csf_filter → ichkq2`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ja` | `bra_csf_index` | 形参 |
| `jb` | `ket_csf_index` | 形参 |
| `k` | `differing_subshell_count` | 变量 |
| `i` | `orbital_index` | 变量 |
| `i_iqa` | `bra_occupation` | 变量 |
| `i_iqb` | `ket_occupation` | 变量 |

保留 `ichkq2`。

## [subshell_occupation_class.f90](../../src/lib/core/subshell_occupation_class.f90)

### `subshell_occupation_class → ichop`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `isubsh` | `subshell_index` | 形参 |
| `icsf` | `csf_index` | 形参 |
| `iocc` | `occupation` | 变量 |
| `ifull` | `subshell_capacity` | 变量 |

保留 `ichop`。

## [integer_vector_copy.f90](../../src/lib/core/integer_vector_copy.f90)

### `integer_vector_copy → icopy`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `n` | `element_count` | 形参 |
| `incx` | `source_stride` | 形参 |
| `incy` | `destination_stride` | 形参 |
| `ix` | `source_values` | 形参 |
| `iy` | `destination_values` | 形参 |
| `m` | `prefix_length` | 变量 |
| `mp1` | `first_full_block_index` | 变量 |

保留 `i`（仅声明）。

## [sparse_eigensolver_initial_basis.f90](../../src/lib/core/sparse_eigensolver_initial_basis.f90)

### `sparse_eigensolver_initial_basis → iniest2`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `nmax` | `maximum_submatrix_order` | 形参 |
| `ncf` | `csf_count` | 形参 |
| `niv` | `requested_eigenvector_count` | 形参 |
| `basis` | `start_vectors_and_eigenvalues` | 形参 |
| `hmx` | `sparse_matrix_values` | 形参 |
| `jcol` | `column_end_offsets` | 形参 |
| `irow` | `row_indices` | 形参 |
| `ns` | `submatrix_order` | 变量 |
| `joffspar` | `previous_sparse_column_end` | 变量 |
| `j` | `column_index` | 变量 |
| `joffnorm` | `packed_column_offset` | 变量 |
| `ir` | `sparse_entry_index` | 变量 |
| `nfound` | `found_eigenvector_count` | 变量 |
| `info` | `lapack_status` | 变量 |
| `ierr` | `negative_lapack_status` | 变量 |
| `iwork` | `integer_workspace` | 变量 |
| `ifail` | `failed_eigenvector_indices` | 变量 |
| `ap` | `packed_submatrix` | 变量 |
| `eigval` | `eigenvalues` | 变量 |
| `vec` | `eigenvectors` | 变量 |
| `work` | `real_workspace` | 变量 |

## [aitken_interpolation.f90](../../src/lib/core/aitken_interpolation.f90)

### `aitken_interpolation → interp`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `narr` | `table_size` | 形参 |
| `xval` | `target_x` | 形参 |
| `yval` | `interpolated_y` | 形参 |
| `accy` | `relative_tolerance` | 形参 |
| `xarr` | `tabulated_x` | 形参 |
| `yarr` | `tabulated_y` | 形参 |
| `mxord` | `maximum_stencil_size` | 常量 |
| `nrstlo` | `lower_bracket_index` | 变量 |
| `nrsthi` | `upper_bracket_index` | 变量 |
| `k` | `search_or_table_index` | 变量 |
| `llo` | `search_lower_bound` | 变量 |
| `lhi` | `search_upper_bound` | 变量 |
| `llr` | `stencil_window_offset` | 变量 |
| `irow` | `stencil_row` | 变量 |
| `locnxt` | `nearest_unused_index` | 变量 |
| `ilirok` | `current_cell` | 变量 |
| `ildiag` | `diagonal_cell` | 变量 |
| `ilothr` | `previous_cell` | 变量 |
| `ibest` | `best_estimate_index` | 变量 |
| `norder` | `stencil_size` | 变量 |
| `dx` | `node_offsets` | 变量 |
| `x` | `stencil_nodes` | 变量 |
| `est` | `estimates` | 变量 |
| `poly` | `aitken_table` | 变量 |
| `diff` | `nearest_offset` | 变量 |
| `difft` | `candidate_offset` | 变量 |
| `debeb` | `best_relative_change` | 变量 |
| `debe` | `estimate_relative_change` | 变量 |
| `used` | `selected_nodes` | 变量 |

### `aitken_interpolation → interp → relative_change`

保留 `change`、`current`、`previous`。

### `aitken_interpolation → interp → iloc`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ind1` | `row_index` | 形参 |
| `ind2` | `column_index` | 形参 |

保留 `iloc`。

## [radial_orbital_interpolation.f90](../../src/lib/core/radial_orbital_interpolation.f90)

### `radial_orbital_interpolation → intrpq`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ma` | `source_grid_size` | 形参 |
| `j` | `orbital_index` | 形参 |
| `dnorm` | `original_norm_squared` | 形参 |
| `pa` | `source_large_component` | 形参 |
| `qa` | `source_small_component` | 形参 |
| `ra` | `source_radial_grid` | 形参 |
| `mxord` | `maximum_stencil_size` | 常量 |
| `i` | `radial_index` | 变量 |
| `mfj` | `last_target_index` | 变量 |
| `nrstlo` | `lower_bracket_index` | 变量 |
| `kount` | `unconverged_point_count` | 变量 |
| `irow` | `stencil_row` | 变量 |
| `k` | `search_or_table_index` | 变量 |
| `nrsthi` | `upper_bracket_index` | 变量 |
| `llo` | `search_lower_bound` | 变量 |
| `lhi` | `search_upper_bound` | 变量 |
| `locnxt` | `nearest_unused_index` | 变量 |
| `ilirok` | `current_cell` | 变量 |
| `ildiag` | `diagonal_cell` | 变量 |
| `ilothr` | `previous_cell` | 变量 |
| `mfjp1` | `first_tail_index` | 变量 |
| `rama` | `source_maximum_radius` | 变量 |
| `x` | `stencil_nodes` | 变量 |
| `dx` | `node_offsets` | 变量 |
| `polyp` | `large_component_table` | 变量 |
| `polyq` | `small_component_table` | 变量 |
| `rn` | `target_maximum_radius` | 变量 |
| `xbar` | `target_radius` | 变量 |
| `pestl` | `previous_large_estimate` | 变量 |
| `qestl` | `previous_small_estimate` | 变量 |
| `diff` | `nearest_offset` | 变量 |
| `difft` | `candidate_offset` | 变量 |
| `dxkmn1` | `previous_node_offset` | 变量 |
| `dxirow` | `current_node_offset` | 变量 |
| `factor` | `inverse_node_separation` | 变量 |
| `pestt` | `large_estimate` | 变量 |
| `qestt` | `small_estimate` | 变量 |
| `dpbp` | `large_relative_change` | 变量 |
| `dqbq` | `small_relative_change` | 变量 |
| `dnfac` | `normalization_factor` | 变量 |
| `set` | `found_candidate` | 变量 |
| `used` | `selected_nodes` | 变量 |

### `radial_orbital_interpolation → intrpq → iloc`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ind1` | `row_index` | 形参 |
| `ind2` | `column_index` | 形参 |

保留 `iloc`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `intrpq` | `def_c::accy` | `accy` | `working_tolerance` |

## [level_selection_parser.f90](../../src/lib/core/level_selection_parser.f90)

### `level_selection_parser → items`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ncmin` | `selected_level_count` | 形参 |
| `ncf` | `maximum_level_index` | 形参 |
| `ierr` | `parse_status` | 形参 |
| `record` | `input_line` | 形参 |
| `ncmin_in` | `initial_level_count` | 变量 |
| `ncd` | `allocated_level_capacity` | 变量 |
| `ifirst` | `range_start_position` | 变量 |
| `istart` | `token_start` | 变量 |
| `i` | `character_index` | 变量 |
| `iend` | `token_end` | 变量 |
| `isize` | `token_width` | 变量 |
| `lenth` | `width_text_length` | 变量 |
| `ios` | `read_status` | 变量 |
| `level` | `parsed_level_index` | 变量 |
| `level1` | `range_start_level` | 变量 |
| `number` | `additional_level_count` | 变量 |
| `j` | `range_offset` | 变量 |
| `reci` | `current_character` | 变量 |
| `form` | `integer_format` | 变量 |
| `cnum` | `width_text` | 变量 |

保留 `length_of_record`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `items` | `def_c::iccmin` | `iccmin` | `selected_level_indices` |

## [angular_momentum_triangle.f90](../../src/lib/core/angular_momentum_triangle.f90)

### `angular_momentum_triangle → itrig`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `i1` | `first_two_j_plus_one` | 形参 |
| `i2` | `second_two_j_plus_one` | 形参 |
| `i3` | `third_two_j_plus_one` | 形参 |
| `i4` | `second_minus_third` | 变量 |

保留 `itrig`。

## [decimal_digit_check.f90](../../src/lib/core/decimal_digit_check.f90)

### `decimal_digit_check → ldigit`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `cst` | `character_to_test` | 形参 |
| `i` | `digit_value` | 变量 |
| `cdgt` | `decimal_digits` | 变量 |

保留 `ldigit`。

## [isotope_data_records.f90](../../src/lib/core/isotope_data_records.f90)

### `isotope_data_records → lodiso`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `a` | `nuclear_mass_number` | 变量 |
| `aparm` | `fermi_diffuseness_fm` | 变量 |
| `cparm` | `fermi_half_density_radius_fm` | 变量 |
| `emnamu` | `nuclear_mass_u` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `lodiso` | `def_c::emn` | `emn` | `nuclear_mass_in_electron_masses` |
| `lodiso` | `nsmdat_c::sqn` | `sqn` | `nuclear_spin` |
| `lodiso` | `nsmdat_c::dmomnm` | `dmomnm` | `nuclear_magnetic_dipole_moment` |
| `lodiso` | `nsmdat_c::qmomb` | `qmomb` | `nuclear_electric_quadrupole_moment` |

## [radial_wavefunction_records.f90](../../src/lib/core/radial_wavefunction_records.f90)

### `radial_wavefunction_records → lodrwf`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `iierr` | `load_status` | 形参 |
| `j` | `orbital_index` | 变量 |
| `i` | `radial_index` | 变量 |
| `k` | `abs_orbital_kappa` | 变量 |
| `nwin` | `loaded_orbital_count` | 变量 |
| `ios` | `read_status` | 变量 |
| `npy` | `file_principal_n` | 变量 |
| `naky` | `file_kappa` | 变量 |
| `my` | `file_grid_size` | 变量 |
| `con` | `charge_times_alpha_squared` | 变量 |
| `fkk` | `kappa_squared` | 变量 |
| `ey` | `file_energy_parameter` | 变量 |
| `pzy` | `file_origin_coefficient` | 变量 |
| `dnorm` | `original_norm_squared` | 变量 |
| `pa` | `file_large_component` | 变量 |
| `qa` | `file_small_component` | 变量 |
| `ra` | `file_radial_grid` | 变量 |

保留 `ierr`（仅声明）。

## [block_level_selection.f90](../../src/lib/core/block_level_selection.f90)

### `block_level_selection → lodstate`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `idblk` | `block_labels` | 形参 |
| `j` | `display_block_index` | 变量 |
| `i` | `selection_index` | 变量 |
| `jblock` | `block_index` | 变量 |
| `ncf` | `block_csf_count` | 变量 |
| `ncminold` | `previous_selection_count` | 变量 |
| `ierr` | `parse_status` | 变量 |
| `ntmp` | `highest_selected_level` | 变量 |
| `str` | `selection_input` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `lodstate` | `def_c::ncmin` | `ncmin` | `selected_level_count` |
| `lodstate` | `def_c::iccmin` | `iccmin` | `selected_level_indices` |
| `lodstate` | `hblock_c::ncmaxblk` | `ncmaxblk` | `highest_selected_level_by_block` |

## [subshell_parent_table.f90](../../src/lib/core/subshell_parent_table.f90)

### `subshell_parent_table → ltab`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `is` | `orbital_indices` | 形参 |
| `nqs` | `occupations` | 形参 |
| `ks` | `subshell_capacities` | 形参 |
| `irows` | `parent_table_rows` | 形参 |
| `kq` | `table_population_index` | 变量 |
| `kq1` | `parent_electron_count` | 变量 |
| `kq2` | `parent_hole_count` | 变量 |
| `i` | `active_shell_position` | 变量 |

## [nuclear_polarization_initialization.f90](../../src/lib/core/nuclear_polarization_initialization.f90)

### `nuclear_polarization_initialization → npol_init`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `path` | `config_path` | 形参 |
| `z` | `nuclear_charge` | 形参 |
| `n` | `grid_point_count` | 形参 |
| `r` | `radial_grid` | 形参 |
| `zz` | `negative_radius_times_potential` | 形参 |
| `cfg` | `config` | 变量 |
| `a` | `dipole_potential` | 变量 |
| `b` | `quadrupole_potential` | 变量 |
| `v` | `polarization_potential` | 变量 |
| `candidate` | `composed_potential` | 变量 |

保留 `status`、`message`。

### `nuclear_polarization_initialization → npol_summary_path`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `slash` | `last_directory_separator` | 变量 |

保留 `path`、`summary`。

### `nuclear_polarization_initialization → npol_write_summary`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `cfg` | `config` | 形参 |
| `n` | `grid_point_count` | 形参 |
| `unit` | `summary_unit` | 变量 |
| `ios` | `io_status` | 变量 |
| `i` | `record_index` | 变量 |
| `bytes` | `file_size_bytes` | 变量 |
| `detail` | `io_message` | 变量 |
| `records` | `summary_records` | 变量 |
| `expected` | `expected_contents` | 变量 |
| `observed` | `read_back_contents` | 变量 |
| `used` | `contents_length` | 变量 |
| `width` | `record_length` | 变量 |

保留 `path`、`status`、`message`、`close_status`。

## [nuclear_polarization_input.f90](../../src/lib/core/nuclear_polarization_input.f90)

### `nuclear_polarization_input → npol_read_config`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `path` | `config_path` | 形参 |
| `cfg` | `config` | 形参 |
| `unit` | `input_unit` | 变量 |
| `ios` | `io_status`、`source_index` | 变量 |
| `bytes` | `input_size_bytes` | 变量 |
| `exists` | `file_exists` | 变量 |
| `text` | `input_contents` | 变量 |
| `detail` | `io_message` | 变量 |

保留 `status`、`message`、`schema_version`、`isotope_z`、`isotope_a`、`enabled`、`use_e1`、`use_e2`、`k1`、`b1`、`k2`、`b2`、`units`、`model_id`、`parameter_source`。

### `nuclear_polarization_input → check_envelope`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `text` | `input_contents` | 形参 |
| `clean` | `normalized_contents` | 变量 |
| `quote` | `open_quote` | 变量 |
| `c` | `current_character` | 变量 |
| `i` | `character_index` | 变量 |
| `last` | `contents_length` | 变量 |
| `finish` | `group_end_index` | 变量 |
| `code` | `character_code` | 变量 |

保留 `status`、`message`。

### `nuclear_polarization_input → npol_validate_config`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `cfg` | `config` | 形参 |
| `z` | `nuclear_charge` | 形参 |

保留 `status`、`message`。

## [nuclear_polarization_potential.f90](../../src/lib/core/nuclear_polarization_potential.f90)

### `nuclear_polarization_potential → npol_check_grid`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `n` | `grid_point_count` | 形参 |
| `r` | `radial_grid` | 形参 |

保留 `status`、`message`。

### `nuclear_polarization_potential → npol_build_grid`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `cfg` | `config` | 形参 |
| `n` | `grid_point_count` | 形参 |
| `r` | `radial_grid` | 形参 |
| `e1` | `dipole_potential` | 形参 |
| `e2` | `quadrupole_potential` | 形参 |
| `total` | `polarization_potential` | 形参 |
| `a` | `dipole_candidate` | 变量 |
| `b` | `quadrupole_candidate` | 变量 |
| `v` | `total_candidate` | 变量 |
| `i` | `radial_index` | 变量 |

保留 `status`、`message`。

### `nuclear_polarization_potential → channel`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `r` | `radius` | 形参 |
| `k` | `channel_strength` | 形参 |
| `b` | `cutoff_radius` | 形参 |
| `p` | `denominator_power` | 形参 |
| `v` | `channel_potential` | 形参 |
| `s` | `radius_scale` | 变量 |
| `m` | `potential_mantissa` | 变量 |
| `e` | `potential_exponent` | 变量 |

保留 `status`、`denominator`。

### `nuclear_polarization_potential → npol_compose`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `n` | `grid_point_count` | 形参 |
| `r` | `radial_grid` | 形参 |
| `coul` | `coulomb_grid` | 形参 |
| `v` | `polarization_potential` | 形参 |
| `total` | `composed_grid` | 形参 |
| `candidate` | `candidate_grid` | 变量 |
| `shift` | `weighted_potential_shift` | 变量 |
| `i` | `radial_index` | 变量 |

保留 `status`、`message`。

## [nuclear_coulomb_potential.f90](../../src/lib/core/nuclear_coulomb_potential.f90)

### `nuclear_coulomb_potential → nucpot`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `i` | `radial_index` | 变量 |
| `nb3` | `points_per_column` | 变量 |
| `nrows` | `table_row_count` | 变量 |
| `ii` | `table_row` | 变量 |
| `ii1` | `first_column_index` | 变量 |
| `ii2` | `second_column_index` | 变量 |
| `ii3` | `third_column_index` | 变量 |
| `c` | `half_density_radius` | 变量 |
| `a` | `diffuseness` | 变量 |
| `abc` | `diffuseness_over_radius` | 变量 |
| `tabc` | `twice_a_over_c` | 变量 |
| `abc2` | `a_over_c_squared` | 变量 |
| `thabc2` | `three_a_over_c_squared` | 变量 |
| `abc3` | `a_over_c_cubed` | 变量 |
| `cba` | `radius_over_diffuseness` | 变量 |
| `pi2` | `pi_squared` | 变量 |
| `hpiac2` | `half_pi_squared_a_over_c_squared` | 变量 |
| `h3php` | `interior_linear_factor` | 变量 |
| `s2mcba` | `origin_inverse_square_sum` | 变量 |
| `s3mcba` | `origin_inverse_cube_sum` | 变量 |
| `sabc3` | `six_a_over_c_cubed` | 变量 |
| `dmsas` | `origin_series_correction` | 变量 |
| `en` | `fermi_normalization` | 变量 |
| `zbn` | `normalized_nuclear_charge` | 变量 |
| `ri` | `radius` | 变量 |
| `rmc` | `radius_from_half_density` | 变量 |
| `rmcba` | `scaled_radius_offset` | 变量 |
| `rbc` | `radius_over_half_density` | 变量 |
| `s2rcba` | `radius_inverse_square_sum` | 变量 |
| `s3rcba` | `radius_inverse_cube_sum` | 变量 |
| `set` | `found_first_exterior_point` | 变量 |

## [file_open.f90](../../src/lib/core/file_open.f90)

### `file_open → openfl`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `nfile` | `file_unit` | 形参 |
| `ierr` | `open_status` | 形参 |
| `filnam` | `file_path` | 形参 |
| `rform` | `record_format` | 形参 |
| `rstat` | `requested_status` | 形参 |
| `ios` | `io_status` | 变量 |
| `loc` | `path_length` | 变量 |

## [radial_orbital_orthogonalization.f90](../../src/lib/core/radial_orbital_orthogonalization.f90)

### `radial_orbital_orthogonalization → orthsc`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `j` | `projected_orbital_indices` | 变量 |
| `l` | `orbital_index` | 变量 |
| `nakl` | `orbital_kappa_value` | 变量 |
| `kount` | `projection_count` | 变量 |
| `mtp0` | `combined_last_grid_index` | 变量 |
| `k` | `previous_orbital_index` | 变量 |
| `mtp` | `last_radial_index` | 变量 |
| `i` | `projection_index` | 变量 |
| `ovlap` | `projection_overlaps` | 变量 |
| `eps` | `tail_cutoff` | 变量 |
| `ovrlap` | `orbital_overlap` | 变量 |
| `dnorm` | `norm_squared` | 变量 |
| `factor` | `normalization_factor` | 变量 |
| `changed` | `projected` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `orthsc` | `def_c::accy` | `accy` | `working_tolerance` |

## [csf_byte_packing.f90](../../src/lib/core/csf_byte_packing.f90)

### `csf_byte_packing → pack_csf_subshell_byte`

`unpacked_value`、`byte_index` 和 `packed_values` 分别表示待写入的值、字节下标和目标数组。现行关键字由[过程指针回归](../../test/f2018/library_modules.f90)覆盖。

## [file_positioning.f90](../../src/lib/core/file_positioning.f90)

### `file_positioning → posfile`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `mode` | `position_mode` | 形参 |
| `nunit` | `file_unit` | 形参 |
| `nrec` | `record_count` | 形参 |
| `forwd` | `rewind_and_advance` | 常量 |
| `backwd` | `backspace_records` | 常量 |
| `i` | `record_index` | 变量 |

## [radial_quadrature.f90](../../src/lib/core/radial_quadrature.f90)

### `radial_quadrature → quad`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `result` | `integral` | 形参 |
| `mtpm1` | `penultimate_grid_index` | 变量 |
| `i` | `radial_index` | 变量 |
| `ip1` | `next_radial_index` | 变量 |
| `tai` | `weighted_value` | 变量 |
| `taip1` | `next_weighted_value` | 变量 |
| `quott` | `weighted_value_ratio` | 变量 |
| `frip1` | `next_integrand` | 变量 |
| `fri` | `integrand` | 变量 |
| `ratio` | `integrand_ratio` | 变量 |
| `rip1` | `next_radius` | 变量 |
| `ri` | `radius` | 变量 |
| `sigma` | `origin_power` | 变量 |

保留 `loc`（仅声明）。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `quad` | `ncc_c::c1` | `c1` | `endpoint_weight` |
| `quad` | `ncc_c::c2` | `c2` | `odd_interior_weight` |
| `quad` | `ncc_c::c3` | `c3` | `middle_weight` |
| `quad` | `ncc_c::c4` | `c4` | `shared_endpoint_weight` |
| `quad` | `tatb_c::ta` | `ta` | `weighted_integrand` |
| `quad` | `tatb_c::mtp` | `mtp` | `last_radial_index` |

## [radial_grid_setup.f90](../../src/lib/core/radial_grid_setup.f90)

### `radial_grid_setup → radgrd`

保留 `grid_points_with_padding`、`radial_index`、`half_grid_point_count`、`debug_row_count`、`debug_row_index`、`left_grid_index`、`right_grid_index`、`exp_grid_step`、`exp_grid_coordinate`、`exp_grid_coordinate_minus_one`、`newton_relative_tolerance`、`grid_step_ratio`、`previous_radius_estimate`、`radius_estimate`、`transformed_coordinate`、`shifted_radius_estimate`、`grid_equation_residual`、`inverse_grid_equation_derivative`、`radius_correction`。

## [weighted_radial_overlap.f90](../../src/lib/core/weighted_radial_overlap.f90)

### `weighted_radial_overlap → rint`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `i` | `first_orbital` | 形参 |
| `j` | `second_orbital` | 形参 |
| `k` | `radial_power` | 形参 |
| `l` | `radial_index` | 变量 |
| `result` | `overlap_integral` | 变量 |

保留 `rint`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `rint` | `tatb_c::ta` | `ta` | `weighted_integrand` |
| `rint` | `tatb_c::mtp` | `mtp` | `last_radial_index` |

## [one_electron_radial_integral.f90](../../src/lib/core/one_electron_radial_integral.f90)

### `one_electron_radial_integral → rinti`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `j` | `bra_orbital` | 形参 |
| `k` | `ket_orbital` | 形参 |
| `mode` | `potential_mode` | 形参 |
| `i` | `radial_index` | 变量 |
| `piece1` | `derivative_term` | 变量 |
| `piece2` | `small_component_mass_term` | 变量 |
| `piece3` | `kappa_term` | 变量 |
| `piece4` | `nuclear_potential_term` | 变量 |

保留 `rinti`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `rinti` | `tatb_c::ta` | `ta` | `integrand_or_large_derivative` |
| `rinti` | `tatb_c::tb` | `tb` | `small_derivative_workspace` |
| `rinti` | `tatb_c::mtp` | `mtp` | `last_radial_index` |

## [rmcdhf_csf_input.f90](../../src/lib/core/rmcdhf_csf_input.f90)

### `rmcdhf_csf_input → rmcdhf_csf_symmetry`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `code` | `packed_symmetry` | 变量 |
| `b` | `block_index` | 变量 |
| `offset` | `cumulative_csf_count` | 变量 |

保留 `csf_index`。

### `rmcdhf_csf_input → descriptor_path`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `name` | `state_name` | 形参 |
| `stem` | `trimmed_name` | 变量 |
| `n` | `name_length` | 变量 |

保留 `path`。

### `rmcdhf_csf_input → load_rmcdhf_descriptors`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `name` | `state_name` | 形参 |
| `ncore` | `core_subshell_count` | 形参 |
| `labels` | `block_labels` | 形参 |
| `source` | `basis_reader` | 变量 |
| `b` | `block_index` | 变量 |
| `first` | `first_csf` | 变量 |
| `last` | `last_csf` | 变量 |
| `offset` | `csf_offset` | 变量 |
| `candidate` | `new_occupations` | 变量 |

保留 `status`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `rmcdhf_descriptor_input` | `syma_c::jpgg` | `jpgg` | `packed_block_symmetries` |

## [physical_constants.f90](../../src/lib/core/physical_constants.f90)

### `physical_constants → setcon`

保留 `reference_bohr_radius_cm`、`inverse_fine_structure_constant`、`reference_speed_of_light_cm_per_s`、`electron_charge_esu`、`electron_mass_g`、`reference_electron_mass_u`、`reference_proton_mass_u`、`reduced_planck_constant_erg_s`、`rydberg_energy_ev`、`rydberg_wavenumber_cm`。

## [csf_basis_loader.f90](../../src/lib/core/csf_basis_loader.f90)

### `csf_basis_loader → setcsla`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `name` | `state_name` | 形参 |
| `ncore` | `core_subshell_count` | 形参 |

保留 `status`。

## [isotope_data_input.f90](../../src/lib/core/isotope_data_input.f90)

### `isotope_data_input → setiso`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `fname` | `isotope_file_path` | 形参 |
| `status` | `requested_status` | 常量 |
| `myname` | `routine_name` | 常量 |
| `form` | `record_format` | 常量 |
| `signature` | `expected_file_header` | 常量 |
| `lenf` | `path_length` | 变量 |
| `ierr` | `open_status` | 变量 |
| `ios` | `read_status` | 变量 |
| `found` | `file_exists` | 变量 |
| `str` | `file_header` | 变量 |

## [angular_recoupling_tables.f90](../../src/lib/core/angular_recoupling_tables.f90)

### `angular_recoupling_tables → setj`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ns` | `active_shell_count` | 形参 |
| `kj23` | `coupling_tables_ready` | 形参 |
| `is` | `orbital_indices` | 形参 |
| `js` | `coupling_node_indices` | 形参 |
| `ks` | `active_two_j_plus_one` | 形参 |
| `ii` | `angular_list_end` | 变量 |
| `ij` | `orbital_index` | 变量 |
| `i` | `table_index` | 变量 |
| `ns1` | `intermediate_coupling_count` | 变量 |
| `ns2` | `coupling_table_size` | 变量 |
| `ns3` | `three_times_shell_count` | 变量 |
| `ij1` | `first_coupling_node` | 变量 |
| `ij2` | `second_coupling_node` | 变量 |

保留 `j`（仅声明）。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `setj` | `m_c::jjc1` | `jjc1` | `bra_intermediate_couplings` |
| `setj` | `m_c::jjc2` | `jjc2` | `ket_intermediate_couplings` |
| `setj` | `m_c::jjq1` | `jjq1` | `bra_subshell_quantum_codes` |
| `setj` | `m_c::jjq2` | `jjq2` | `ket_subshell_quantum_codes` |
| `setj` | `m_c::jlist` | `jlist` | `peel_orbital_indices` |

## [machine_constants.f90](../../src/lib/core/machine_constants.f90)

### `machine_constants → setmc`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `dnum` | `machine_limit` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `setmc` | `def_c::expmax` | `expmax` | `maximum_exp_argument` |
| `setmc` | `def_c::expmin` | `expmin` | `minimum_exp_argument` |
| `setmc` | `def_c::precis` | `precis` | `machine_precision` |
| `setmc` | `def_c::tenmax` | `tenmax` | `maximum_decimal_exponent` |

## [radial_dirac_coefficients.f90](../../src/lib/core/radial_dirac_coefficients.f90)

### `radial_dirac_coefficients → setpot`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `j` | `orbital_index` | 形参 |
| `jp` | `matching_grid_index` | 形参 |
| `i` | `radial_index` | 变量 |
| `dmhb2c` | `negative_step_over_twice_c` | 变量 |
| `energy` | `binding_energy` | 变量 |
| `enefac` | `twice_c_squared_minus_energy` | 变量 |
| `rpi` | `grid_derivative` | 变量 |
| `yprpor` | `scaled_orbital_potential` | 变量 |
| `notset` | `searching_for_match` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `setpot` | `def_c::accy` | `accy` | `working_tolerance` |
| `setpot` | `int_c::tf` | `tf` | `large_equation_coefficients` |
| `setpot` | `int_c::tg` | `tg` | `small_equation_coefficients` |
| `setpot` | `pote_c::yp` | `yp` | `orbital_y_potential` |

## [radial_stencil_weights.f90](../../src/lib/core/radial_stencil_weights.f90)

### `radial_stencil_weights → setqic`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `i` | `weight_index` | 变量 |
| `j` | `column_index` | 变量 |
| `b13` | `derivative_numerators` | 变量 |
| `cg` | `starting_weights` | 变量 |
| `c5num` | `five_point_numerators` | 变量 |
| `c6num` | `six_point_numerators` | 变量 |
| `b13den` | `derivative_denominator` | 变量 |
| `denom` | `starting_denominator` | 变量 |
| `c5den` | `five_point_denominator` | 变量 |
| `c6den` | `six_point_denominator` | 变量 |
| `factor` | `weight_scale` | 变量 |
| `first` | `initialize_derivative_table` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `setqic` | `cnc_c::cnc5c` | `cnc5c` | `five_point_weights` |
| `setqic` | `cnc_c::cnc6c` | `cnc6c` | `six_point_weights` |
| `setqic` | `lic13_c::a` | `a13` | `derivative_weights` |
| `setqic` | `ncc_c::c1` | `c1` | `endpoint_weight` |
| `setqic` | `ncc_c::c2` | `c2` | `odd_interior_weight` |
| `setqic` | `ncc_c::c3` | `c3` | `middle_weight` |
| `setqic` | `ncc_c::c4` | `c4` | `shared_endpoint_weight` |
| `setqic` | `sbc_c::c` | `c` | `starting_step_weights` |

## [csf_pair_quantum_numbers.f90](../../src/lib/core/csf_pair_quantum_numbers.f90)

### `csf_pair_quantum_numbers → setqna`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ja` | `bra_csf_index` | 形参 |
| `jb` | `ket_csf_index` | 形参 |
| `j` | `orbital_or_peel_index` | 变量 |
| `k` | `quantum_number_selector` | 变量 |
| `jcnt` | `next_coupling_index` | 变量 |
| `jcntop` | `open_shell_count` | 变量 |
| `jw1` | `first_peel_orbital` | 变量 |
| `jw2` | `second_peel_orbital` | 变量 |
| `jw` | `peel_orbital` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `setqna` | `m_c::jjc1` | `jjc1` | `bra_intermediate_couplings` |
| `setqna` | `m_c::jjc2` | `jjc2` | `ket_intermediate_couplings` |
| `setqna` | `m_c::jjq1` | `jjq1` | `bra_subshell_quantum_codes` |
| `setqna` | `m_c::jjq2` | `jjq2` | `ket_subshell_quantum_codes` |
| `setqna` | `m_c::jlist` | `jlist` | `peel_orbital_indices` |
| `setqna` | `m_c::klist` | `klist` | `core_orbital_indices` |
| `setqna` | `m_c::ncore` | `ncore` | `core_orbital_count` |
| `setqna` | `m_c::npeel` | `npeel` | `peel_orbital_count` |
| `setqna` | `m_c::nq1` | `nq1` | `bra_occupations` |
| `setqna` | `m_c::nq2` | `nq2` | `ket_occupations` |

## [radial_wavefunction_input.f90](../../src/lib/core/radial_wavefunction_input.f90)

### `radial_wavefunction_input → setrwfa`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `name` | `radial_file_path` | 形参 |
| `ierr` | `load_status` | 变量 |
| `ios` | `read_status` | 变量 |
| `g92rwf` | `file_signature` | 变量 |

## [coulomb_tensor_rank_range.f90](../../src/lib/core/coulomb_tensor_rank_range.f90)

### `coulomb_tensor_rank_range → skrc`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `kd1` | `first_direct_rank` | 形参 |
| `kd2` | `direct_rank_count` | 形参 |
| `ke1` | `first_exchange_rank` | 形参 |
| `ke2` | `exchange_rank_count` | 形参 |
| `is` | `orbital_indices` | 形参 |
| `kaps` | `orbital_kappas` | 形参 |
| `ks` | `orbital_two_j_plus_one` | 形参 |
| `isd1` | `direct_kappa_sign_a` | 变量 |
| `isd2` | `direct_kappa_sign_b` | 变量 |
| `kd1a` | `twice_direct_lower_a` | 变量 |
| `kd1b` | `twice_direct_lower_b` | 变量 |
| `kd2a` | `twice_direct_upper_a` | 变量 |
| `kd2b` | `twice_direct_upper_b` | 变量 |
| `ise1` | `exchange_kappa_sign_a` | 变量 |
| `ise2` | `exchange_kappa_sign_b` | 变量 |
| `ke1a` | `twice_exchange_lower_a` | 变量 |
| `ke1b` | `twice_exchange_lower_b` | 变量 |
| `ke2a` | `twice_exchange_upper_a` | 变量 |
| `ke2b` | `twice_exchange_upper_b` | 变量 |

## [slater_radial_integral.f90](../../src/lib/core/slater_radial_integral.f90)

### `slater_radial_integral → slater`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ia` | `orbital_a` | 形参 |
| `ib` | `orbital_b` | 形参 |
| `ic` | `orbital_c` | 形参 |
| `id` | `orbital_d` | 形参 |
| `k` | `tensor_rank` | 形参 |
| `i` | `radial_index` | 变量 |
| `result` | `slater_integral` | 变量 |

保留 `slater`。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `slater` | `tatb_c::ta` | `ta` | `weighted_integrand` |
| `slater` | `tatb_c::tb` | `tb` | `hartree_y_values` |
| `slater` | `tatb_c::mtp` | `mtp` | `last_radial_index` |

## [mcp_coefficient_output.f90](../../src/lib/core/mcp_coefficient_output.f90)

### `mcp_coefficient_output → speak`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ja` | `bra_csf_index` | 形参 |
| `jb` | `ket_csf_index` | 形参 |
| `ia1` | `bra_orbital_1` | 形参 |
| `ib1` | `bra_orbital_2` | 形参 |
| `ia2` | `ket_orbital_1` | 形参 |
| `ib2` | `ket_orbital_2` | 形参 |
| `k` | `tensor_rank` | 形参 |
| `x` | `coefficient` | 形参 |

## [sparse_symmetric_matvec.f90](../../src/lib/core/sparse_symmetric_matvec.f90)

### `sparse_symmetric_matvec → spicmv2`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `n` | `matrix_order` | 形参 |
| `m` | `vector_count` | 形参 |
| `b` | `input_vectors` | 形参 |
| `c` | `product_vectors` | 形参 |
| `ibeg` | `column_start` | 变量 |
| `icol` | `column_index` | 变量 |
| `iend` | `column_end` | 变量 |
| `nelc` | `column_entry_count` | 变量 |
| `iv` | `vector_index` | 变量 |
| `diag` | `diagonal_plus_previous_terms` | 变量 |
| `dl` | `off_diagonal_dot_product` | 变量 |

## [radial_origin_solution.f90](../../src/lib/core/radial_origin_solution.f90)

### `radial_origin_solution → start`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `mxiter` | `maximum_iterations` | 常量 |
| `kiorb` | `orbital_kappa_value` | 变量 |
| `i` | `radial_or_coefficient_index` | 变量 |
| `jorb` | `direct_orbital_index` | 变量 |
| `niter` | `iteration_count` | 变量 |
| `j` | `start_grid_index` | 变量 |
| `rdp` | `large_equation_coupling` | 变量 |
| `rdq` | `small_equation_coupling` | 变量 |
| `rsep` | `large_equation_source` | 变量 |
| `rseq` | `small_equation_source` | 变量 |
| `spest` | `large_component_estimates` | 变量 |
| `sqest` | `small_component_estimates` | 变量 |
| `obc` | `inverse_speed_of_light` | 变量 |
| `zonc` | `charge_over_c` | 变量 |
| `giorb` | `origin_exponent` | 变量 |
| `fkiorb` | `kappa_real` | 变量 |
| `omgi` | `one_minus_origin_exponent` | 变量 |
| `omgmk` | `one_minus_exponent_minus_kappa` | 变量 |
| `omgpk` | `one_minus_exponent_plus_kappa` | 变量 |
| `rsep1` | `large_origin_source_sum` | 变量 |
| `rseq1` | `small_origin_source_sum` | 变量 |
| `pzero` | `direct_orbital_origin_parameter` | 变量 |
| `p1` | `large_origin_value` | 变量 |
| `q1` | `small_origin_value` | 变量 |
| `sump` | `large_equation_sum` | 变量 |
| `sumq` | `small_equation_sum` | 变量 |
| `csq` | `c_squared` | 变量 |
| `twocsq` | `twice_c_squared` | 变量 |
| `energy` | `binding_energy` | 变量 |
| `enefac` | `twice_c_squared_minus_energy` | 变量 |
| `ri` | `radius` | 变量 |
| `rpi` | `grid_derivative` | 变量 |
| `rirpi` | `radius_times_derivative` | 变量 |
| `ypirpi` | `potential_times_derivative` | 变量 |
| `difmaw` | `unscaled_maximum_change` | 变量 |
| `difmax` | `maximum_change` | 变量 |
| `rj` | `start_grid_radius` | 变量 |

保留 `orbital_index`、`radial_equation_type`、`origin_slope_parameter`、`small_component_origin_coefficient`、`large_component_values`、`small_component_values`、`factor`、`coefij`（仅声明）、`pi`（仅声明）、`qi`（仅声明）。

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `start` | `cnc_c::cnc6c` | `cnc6c` | `six_point_weights` |
| `start` | `def_c::accy` | `accy` | `working_tolerance` |
| `start` | `pote_c::yp` | `yp` | `orbital_y_potential` |
| `start` | `pote_c::xp` | `xp` | `large_exchange_source` |
| `start` | `pote_c::xq` | `xq` | `small_exchange_source` |
| `start` | `scf_c::ndcof` | `ndcof` | `direct_coefficient_count` |
| `start` | `scf_c::nda` | `nda` | `direct_orbital_indices` |
| `start` | `scf_c::da` | `da` | `direct_coefficients` |

## [run_timer_start.f90](../../src/lib/core/run_timer_start.f90)

### `run_timer_start → starttime`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ncount1` | `start_clock_count` | 形参 |
| `progname` | `program_name` | 形参 |
| `ncount_rate` | `clock_tick_rate` | 变量 |
| `ncount_max` | `maximum_clock_count` | 变量 |

保留 `nymduhmsm`（仅声明）、`chdate`（仅声明）、`chtime`（仅声明）、`chzone`（仅声明）、`msg`（仅声明）。

## [run_timer_finish.f90](../../src/lib/core/run_timer_finish.f90)

### `run_timer_finish → stoptime`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `ncount1` | `start_clock_count` | 形参 |
| `progname` | `program_name` | 形参 |
| `ncount2` | `elapsed_clock_count` | 变量 |
| `ncount_rate` | `clock_tick_rate` | 变量 |
| `ncount_max` | `maximum_clock_count` | 变量 |
| `nseconds` | `elapsed_seconds` | 变量 |
| `nymduhmsm` | `date_time_values` | 变量 |
| `chdate` | `date_text` | 变量 |
| `chtime` | `time_text` | 变量 |
| `chzone` | `time_zone_text` | 变量 |
| `str2nds` | `elapsed_seconds_text` | 变量 |
| `msg` | `elapsed_time_message` | 变量 |

## [posix_services.f90](../../src/lib/core/posix_services.f90)

### `posix_services → c_exit`

保留 `status`。

### `posix_services → c_system`

保留 `command`、`status`。

### `posix_services → c_chdir`

保留 `path`、`status`。

### `posix_services → c_getcwd`

保留 `buffer`、`capacity`、`address`。

### `posix_services → c_access`

保留 `path`、`mode`、`status`。

### `posix_services → terminate_process`

保留 `status`。

### `posix_services → run_shell`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `buffer` | `c_command` | 变量 |
| `i` | `character_index` | 变量 |
| `n` | `command_length` | 变量 |

保留 `status`、`command`。

### `posix_services → change_directory`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `buffer` | `c_path` | 变量 |
| `i` | `character_index` | 变量 |
| `n` | `path_length` | 变量 |

保留 `status`、`path`。

### `posix_services → working_directory`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `buffer` | `c_path` | 变量 |
| `address` | `buffer_address` | 变量 |
| `i` | `character_index` | 变量 |

保留 `status`、`path`。

### `posix_services → path_access`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `buffer` | `c_path` | 变量 |
| `mask` | `access_mask` | 变量 |
| `i` | `character_index` | 变量 |
| `n` | `path_length` | 变量 |

保留 `status`、`path`、`mode`。

### `posix_services → make_directory`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `prefix` | `command_prefix` | 常量 |
| `i` | `character_index` | 变量 |
| `used` | `command_length` | 变量 |

保留 `path`、`status`、`command`。

## [hartree_yz_functions.f90](../../src/lib/core/hartree_yz_functions.f90)

### `hartree_yz_functions → yzk`

| 原名 | 迁移后 | 角色 |
| --- | --- | --- |
| `k` | `tensor_rank` | 形参 |
| `i` | `first_orbital` | 形参 |
| `j` | `second_orbital` | 形参 |
| `klimit` | `maximum_tensor_rank` | 常量 |
| `klimit1` | `rank_cache_size` | 常量 |
| `ii` | `radial_index` | 变量 |
| `mtpp1` | `first_padding_index` | 变量 |
| `mtpp3` | `last_quadrature_index` | 变量 |
| `mtpp4` | `first_asymptotic_index` | 变量 |
| `kk` | `stencil_or_tail_index` | 变量 |
| `np4` | `padded_grid_end` | 变量 |
| `rhop` | `weighted_orbital_product` | 变量 |
| `rttk` | `radius_to_rank` | 变量 |
| `rttkm` | `inverse_radius_to_rank` | 变量 |
| `rttk1` | `radius_to_rank_plus_one` | 变量 |
| `rttkm1` | `inverse_radius_to_rank_plus_one` | 变量 |
| `rm` | `inverse_radius` | 变量 |
| `wk` | `hartree_outer_integral` | 变量 |
| `temp` | `weighted_integrand` | 变量 |
| `sum` | `stencil_integral` | 变量 |
| `rtmp` | `interval_integral` | 变量 |
| `zklim` | `inner_integral_limit` | 变量 |
| `dif` | `tail_correction` | 变量 |
| `kcalc` | `rank_cached` | 变量 |

### 共享状态在本文件中的别名

| 作用域 | 定义 | 原局部名称 | 迁移后 |
| --- | --- | --- | --- |
| `yzk` | `cnc_c::cnc5c` | `cnc5c` | `five_point_weights` |
| `yzk` | `def_c::accy` | `accy` | `working_tolerance` |
| `yzk` | `ncc_c::c1` | `c1` | `endpoint_weight` |
| `yzk` | `ncc_c::c2` | `c2` | `odd_interior_weight` |
| `yzk` | `ncc_c::c3` | `c3` | `middle_weight` |
| `yzk` | `tatb_c::ta` | `zk` | `hartree_z_values` |
| `yzk` | `tatb_c::tb` | `yk` | `hartree_y_values` |
| `yzk` | `tatb_c::mtp` | `mtp` | `last_radial_index` |
