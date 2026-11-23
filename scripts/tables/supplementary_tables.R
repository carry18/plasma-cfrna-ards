# 输入：固定Table 1临床CSV、02–06分析结果及FigureS2生成的重采样汇总。
# 输出：Table 1及Supplementary Tables S1–S5。
# 在项目根目录运行；可用CF_PYTHON指定Python解释器。
python <- Sys.getenv("CF_PYTHON", "python3")
status <- system2(python, "scripts/internal/supplementary_tables.py")
if (status != 0L) stop("Table generation failed")
