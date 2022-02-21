> library(MatrixEQTL)
> SNP_file_name = "D:/dc/项目/结直肠腺瘤/RNA部分分析结果/19Samples_snp_matrixV2.txt"
> snps = SlicedData$new();
> snps$fileDelimiter = "\t";
> snps$fileOmitCharacters = "NA";
> snps$fileSkipRows = 1; 
> snps$fileSkipColumns = 1;
> snps$fileSliceSize = 2000;
> snps$LoadFile( SNP_file_name );
> expression_file_name = "D:/dc/项目/结直肠腺瘤/RNA部分分析结果/19samples_crcRNAtpm.txt"
> gene = SlicedData$new();
> gene$fileDelimiter = "\t";
> gene$fileOmitCharacters = "NA";
> gene$fileSkipRows = 1; 
> gene$fileSkipColumns = 1; 
> gene$fileSliceSize = 2000;
> gene$LoadFile(expression_file_name);
Rows read: 504 done.
> covariates_file_name = "D:/dc/项目/结直肠腺瘤/RNA部分分析结果/covariates_file.txt"
> cvrt = SlicedData$new();
> cvrt$fileDelimiter = "\t";
> cvrt$fileOmitCharacters = "NA";
> cvrt$fileSkipRows = 1;
> cvrt$fileSkipColumns = 1;
> cvrt$fileSliceSize = 2000;
> cvrt$LoadFile(covariates_file_name);
> snps_location_file_name = "D:/dc/项目/结直肠腺瘤/RNA部分分析结果/19Samples_snp_coordV1.txt"
> snpspos = read.table(snps_location_file_name,header = TRUE,stringsAsFactors = FALSE);
> gene_location_file_name = "D:/dc/项目/结直肠腺瘤/RNA部分分析结果/gene_position.txt"
> genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
> useModel = modelLINEAR;
> output_file_name_cis = "D:/dc/项目/结直肠腺瘤/RNA部分分析结果/cis-eqtl"
> output_file_name_tra = "D:/dc/项目/结直肠腺瘤/RNA部分分析结果/trans-eqtl"
> pvOutputThreshold_cis = 2e-2;
> pvOutputThreshold_tra = 1e-2;
> errorCovariance = numeric();
> cisDist = 1e6;
> me = Matrix_eQTL_main(
+     snps = snps, 
+     gene = gene, 
+     cvrt = cvrt,
+     output_file_name     = output_file_name_tra,
+     pvOutputThreshold     = pvOutputThreshold_tra,
+     useModel = useModel, 
+     errorCovariance = errorCovariance, 
+     verbose = TRUE, 
+     output_file_name.cis = output_file_name_cis,
+     pvOutputThreshold.cis = pvOutputThreshold_cis,
+     snpspos = snpspos, 
+     genepos = genepos,
+     cisDist = cisDist,
+     pvalue.hist = "qqplot",
+     min.pv.by.genesnp = FALSE,
+     noFDRsaveMemory = FALSE);
> unlink(output_file_name_tra);
> unlink(output_file_name_cis);
> cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
> cat('Detected local eQTLs:', '\n');
> cat('Detected distant eQTLs:', '\n');
> show(me$trans$eqtls)
> show(me$cis$eqtls)
> plot(me)
> write.table(me$cis$eqtls,"D:/dc/项目/结直肠腺瘤/RNA部分分析结果/cis-eqtl.txt",sep="\t")
> write.table(me$trans$eqtls,"D:/dc/项目/结直肠腺瘤/RNA部分分析结果/trans-eqtl.txt",sep="\t")




