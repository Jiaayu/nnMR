# nnMR
## 一些常用的函数
### 486种代谢物批量分析
示例
met_486(exposure_path='D:\\486\\已处理',outcome_path='finngen_R9_DM_NEPHROPATHY_EXMORE',
        p_threshold=5e-6,kb=10000,r2=0.001,outcome_source = 'finn')
exposure_path:暴露所在的路径 文件需要为csv
outcome_path:结局所在的路径
p_threshold,kb,r2为筛选和clump参数
outcome_source可选"finn"和"ieu"

### 1400种代谢物批量分析
met_1400(exposure_path='D:\\1400\\exposure\\已处理',outcome_path='finngen_R9_DM_NEPHROPATHY_EXMORE',
         p_threshold=1e-6,kb=10000,r2=0.001,outcome_source = 'finn')
同上

### LDSC批量分析
nn_LDSC_auto(exposure_path='D:\\1400\\原文件',exposure_source=1400,
             outcome_path='finngen_R9_K11_GASTRODUOULC.gz',outcome_source='finn',outcome_case=329603)
exposure_source可选1400、486或rewrite
outcome_case根据查到的信息进行填写
