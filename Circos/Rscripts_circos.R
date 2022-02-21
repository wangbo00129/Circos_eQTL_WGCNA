library(RCircos)
> RCircos.Set.Plot.Area()
> RCircos.Highligh.Chromosome.Ideogram(highlight.width = 8,highlight.pos=1.8)
> RCircos.Label.Chromosome.Names(chr.name.pos=NULL)
> cyto.info = read.table("D:/dc/项目/结直肠腺瘤/RNA部分分析结果/25Ade_UCSC.HG19.Human.CytoBandIdeogram",sep="\t",header=T,row.names = 1)
> RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,tracks.inside=3, tracks.outside=0 )
> RCircos.Set.Plot.Area()
> RCircos.Chromosome.Ideogram.Plot()
> RCircos.Link.Data = read.table("D:/dc/项目/结直肠腺瘤/RNA部分分析结果/25Ade_RCircos.Link.Data",sep="\t",header=T,row.names = 1)
> RCircos.Link.Plot(RCircos.Link.Data,track.num=1,by.chromosome=FALSE)





