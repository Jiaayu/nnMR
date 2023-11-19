# nnMR
## 安装及载入nnMR
1. 如果没有安装devtools则安装devtools
```
install.packages("devtools")
```
2. 安装nnMR
```
devtools::install_github("Jiaayu/nnMR")
```
3. 载入nnMR
```
library(nnMR)
```
## 一些常用的函数
### 一、486种代谢物批量分析
示例
~~met_486(exposure_path='D:\\486\\已处理',outcome_path='finngen_R9_DM_NEPHROPATHY_EXMORE',
        p_threshold=5e-6,kb=10000,r2=0.001,outcome_source = 'finn')
exposure_path:暴露所在的路径 文件需要为csv
outcome_path:结局所在的路径
p_threshold,kb,r2为筛选和clump参数
outcome_source可选"finn"和"ieu"~~

### 二、1400种代谢物批量分析
~~met_1400(exposure_path='D:\\1400\\exposure\\已处理',outcome_path='finngen_R9_DM_NEPHROPATHY_EXMORE',
         p_threshold=1e-6,kb=10000,r2=0.001,outcome_source = 'finn')~~
同上

### 三、LDSC批量分析
~~nn_LDSC_auto(exposure_path='D:\\1400\\原文件',exposure_source=1400,
             outcome_path='finngen_R9_K11_GASTRODUOULC.gz',outcome_source='finn',outcome_case=329603)
exposure_source可选1400、486或rewrite
outcome_case根据查到的信息进行填写~~


### 四、731种免疫循环批量分析
```
suppressMessages(immnune_auto(exposure_path='D:/immune/csv',outcome_data=outcome_data,output_path='result.csv',
                 multiprocess=FALSE,start_num=NA,end_num=NA,
                 exposure_pval=1e-5,clump_r2=0.1,clump_kb=500,bfile_path='D:/clump_pop/EUR'))
```

1. exposure_path:为731种免疫细胞的csv文件所在目录
2. outcome_data:必须为数据框，为结局数据，**请看以下结局文件的处理**
3. output_path:为输出结果文件的名字，默认为result.csv，如果使用多线程分析则每个线程输出文件名字必须不同（否则会覆盖）
4. multiprocess:是否启用多线程，默认关闭
5. start_num, end_num:多线程的分析起点和终点，如果要使用多线程，则分别设置每个线程的起点和终点，比如线程1设置为1和350，线程2设置为351和731（总共731）
6. exposure_pval:筛选暴露的P值，默认为1e-5
7. clump_r2,clump_kb:clump使用的r2和kb，默认为0.1和500
8. bfile_path:clump所需要的文件路径


### 五、731种免疫循环在线批量分析
```
suppressMessages(immnune_online(outcome_data,output_path='result.csv',
                 multiprocess=FALSE,start_num=NA,end_num=NA,
                 exposure_pval=1e-5,clump_r2=0.1,clump_kb=500))
```

1. outcome_data:必须为数据框格式的结局数据
3. output_path:为输出结果文件的名字，默认为result.csv，如果使用多线程分析则每个线程输出文件名字必须不同（否则会覆盖）
4. multiprocess:是否启用多线程，默认关闭
5. start_num, end_num:多线程的分析起点和终点，如果要使用多线程，则分别设置每个线程的起点和终点，比如线程1设置为1和350，线程2设置为351和731（总共731）
6. exposure_pval:筛选暴露的P值，默认为1e-5
7. clump_r2,clump_kb:clump使用的r2和kb，默认为0.1和500

### 六、731种免疫循环在线 <反向> 批量分析
```
suppressMessages(immnune_online_reverse(exposure_data,output_path='result_reverse.csv',
                                        multiprocess=FALSE,start_num=NA,end_num=NA,
                                        exposure_pval=5e-8,clump_r2=0.001,clump_kb=10000'))
```

1. exposure_data:必须为数据框格式的暴露数据，**注意format函数把type改为exposure**
3. output_path:为输出结果文件的名字，默认为result_reverse.csv，如果使用多线程分析则每个线程输出文件名字必须不同（否则会覆盖）
4. multiprocess:是否启用多线程，默认关闭
5. start_num, end_num:多线程的分析起点和终点，如果要使用多线程，则分别设置每个线程的起点和终点，比如线程1设置为1和350，线程2设置为351和731（总共731）
6. exposure_pval:筛选暴露的P值，默认为5e-8
7. clump_r2,clump_kb:clump使用的r2和kb，默认为0.001和10000

#### 循环结束后进行ID-表型转换
```
transform_immune_process(dat='result.csv')
```

#### 循环结束后进行FDR校正
```
nnFDR(filepath='result.csv')
```


## 结局文件的处理（以芬兰数据库为例）
1. 读取结局
2. 增加一列samplesize，值为样本量
3. 增加一列phenotype，值为结局的疾病名
4. 示例代码：
```
out_dat <- data.table::fread('finngen_R9_L12_ATOPIC.gz')
out_dat$samplesize <- 350062
out_dat$phenotype <- 'Atopic dermatitis'
out_dat <- TwoSampleMR::format_data(dat = out_dat,type = 'outcome',
                                    snp_col = 'rsids',chr_col = '#chrom',pos_col = 'pos',
                                    effect_allele_col = 'alt',other_allele_col = 'ref',
                                    pval_col = 'pval',beta_col = 'beta',se_col = 'sebeta',
                                    eaf_col = 'af_alt',samplesize_col = 'samplesize',phenotype_col = 'phenotype')
```

## 以下是不同GWAS数据库的格式化代码，方便自用
##### 芬兰数据库
```
data <- TwoSampleMR::format_data(dat = dat,type = 'outcome',
                                    snp_col = 'rsids',chr_col = '#chrom',pos_col = 'pos',
                                    effect_allele_col = 'alt',other_allele_col = 'ref',
                                    pval_col = 'pval',beta_col = 'beta',se_col = 'sebeta',
                                    eaf_col = 'af_alt')
```
##### 731免疫细胞
```
data <- TwoSampleMR::format_data(dat = dat,type = 'exposure',
                                   chr_col = 'chromosome',pos_col = 'base_pair_location',
                                   effect_allele_col = 'effect_allele',other_allele_col = 'other_allele',
                                   samplesize_col = 'n',eaf_col = 'effect_allele_frequency',
                                   beta_col = 'beta',se_col = 'standard_error',pval_col = 'p_value',snp_col = 'variant_id')
```
