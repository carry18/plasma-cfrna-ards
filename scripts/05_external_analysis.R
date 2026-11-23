# 输入：prepare生成的source_paths.json、当前网络analysis.rds，以及GEO原始文件。
# 输出：GEO基因表达/元数据、外部评分及成员相关结果、更新后的网络analysis.rds。
# 在项目根目录运行；前置步骤：01、02、03。
python <- Sys.getenv("CF_PYTHON", "python3")
status <- system2(python, "scripts/internal/geo_parse.py")
if (status != 0L) stop("GEO parsing failed")
sys.source("scripts/internal/geo_prepare.R", envir = new.env())
source("scripts/internal/external_analysis.R")
