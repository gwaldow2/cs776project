library(PROPER)
sim.opts.Cheung = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
                                         lOD="cheung", lBaselineExpr="cheung")
sim.opts.Bottomly = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
                                           lOD="bottomly", lBaselineExpr="bottomly")
sim.opts.MAQC = RNAseq.SimOptions.2grp(ngenes = 20000, p.DE=0.05,
                                           lOD="maqc", lBaselineExpr="maqc")


d <- simRNAseq(sim.opts.Bottomly, n1=3, n2=3)


sim.opts.Cheung$lOD
sim.opts.Cheung$lBaselineExpr