# 木星黎明风暴论文代码包

拟发布版本：**v1.1.0**。整理日期：2026-10-18。本地文件已整理，不代表已上传。
按作者选择，暂不指定许可证；没有添加 MIT、CC BY 或其他授权声明。

这份包整合当前 v4 正文和 SM 使用的主要分析代码，限定于 **2017-05-16**
黎明风暴事件和 **2017-05-10** 宁静对照，不含2019年HST项目。

## 最快开始

解压后，进入含 `run_analysis.py` 的文件夹，在地址栏输入 `powershell` 并回车：

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_analysis.py check
.\.venv\Scripts\python.exe run_analysis.py figures
```

图片在 `output/figures/`。命令直接使用包内处理后观测数据，不要求拥有原作者硬盘路径。
不用双击 `.py` 文件，也不用激活虚拟环境。

## 重新计算

```powershell
.\.venv\Scripts\python.exe run_analysis.py psd
.\.venv\Scripts\python.exe run_analysis.py figures
```

`psd` 从三分量磁场重新分解、去均值、计算小波、排除边界影响并标定PSD，同时做数值验证；
不是简单读入已有PSD。可能需要数分钟。随后 `figures` 将使用输出目录中的新PSD绘图。

```powershell
.\.venv\Scripts\python.exe run_analysis.py cwt
.\.venv\Scripts\python.exe run_analysis.py bandpass
.\.venv\Scripts\python.exe run_analysis.py theory
```

分别为独立重算小波振幅、S1带通分析、运行原始理论求根脚本。
`all` 可以运行所有流程。生成文件统一放在 `output/`，随包参考输入不会被覆盖。

## 这次比旧包新增了什么

- 最新单边小波PSD标定、验证代码和完整谱表。
- Fig4、S3中 3.5×10⁻³–5×10⁻² Hz 展示范围对应的绘图代码。
- Fig3和S3原始小波振幅显示数组，以及2017年处理后磁场窗口。
- 可移植入口、依赖表、文件校验、变量单位、代码与图号对应关系。
- GitHub上传、发布Release和Zenodo更新版本的中文步骤。

完整谱仍保留约3.359 mHz的平行分量局部峰；显示区间缩窄不会删除原数据。
“展示频段内平行无局部峰”不等于“全频段平行无峰”或“垂直功率占优”。

## 公开前后

按 [公开步骤_中文.md](docs/公开步骤_中文.md) 操作。根据你此前提供的引用，已有
GitHub/Zenodo v1.0.1；优先更新原仓库并归档新版本。原DOI不直接改写成新版本DOI。
发布后把新版本页面实际生成的DOI更新到正文参考文献[52]和引用文件中。

本包没有声称完成原始PDS数据到JSS/RTP的坐标转换，也不含HST数据处理及极光制图代码。
理论脚本按原样保留，最终Fig1曲线分支的精确对应仍需单独核实。
