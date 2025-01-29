# for convenience, create an object in memory that is only the RNAseq data (no metadata), rows are individuals, columns are the many different genes
myRNAdata <-   mydataset [ ,  columns_w_RNA_data ]

# Create an object in memory that is only the metadata

metadata <-  mydataset [ ,  columns_w_metadata ]

# Build a data table to insert results into as you proceed through a loop
results.table <- as.data.frame(   matrix( nrow = #_genes,  ncol = 10 ) )
                                            
                                            # iterate through all your genes
                                            for( i in 1:ncol(myRNAdata)  {
                                              # for each gene grab its read count vector
                                              focal_gene <- myRNAdata[ , i ]
                                              # get the sum of all other genes (the reads that are NOT the focal gene for this iteration)
                                              all_others <-  rowSums( myRNAdata[ , -i ]
                                                                      # create a dependent variable Y that is 2 columns, one is the counts of the focal gene, the other is the counts of everythng else. This is a binomial output of a list of “successes” and “failures” (counting the focal gene, or not)
                                                                      Y <- cbind(focal_gene, all_others)
                                                                      # Now for the GLM
                                                                      Model <- glm( Y ~  Population * treatment, family = “quasipoisson”)   
                                                                      # I tend to favor quasipoisson over Poisson because it is more conservative, less sensitive to high leverage outliers
                                                                      # Now extract the results and save in the results.table
                                                                      gene_results  <- summary( Model)
                                                                      # or a ANOVA version
                                                                      gene_results_2 <- anova(Model, test = “LRT”)
                                                                      # Now populate the results.table row i with whatever information you want to save about this iteration’s focal gene:
                                                                      results.table[ i , 1]  <- colnames(myRNAdata)[i]    # saves gene ID or name
                                                                      results.table[i , 2]  <- sum(myRNAdata[, i] / sum(myRNAdata)   # saves mean relative expression level, for reference - really rare genes might not be of interest and you can sort them and exclude rare ones later
                                                                                                   results.table[i, 3] <-  gene_results$coef[2,4]  # saves P value for Population effect; you may need to check the details here to suit your model
                                                                                                   results.table[i, 4] <-  gene_results$coef[3,4]  # saves P value for Treatment effect
                                                                                                   results.table[i, 5] <-  gene_results$coef[4,4]  # saves P value for Population*treatment interaction effect
                                                                                                   # you might save additional information: means for each treatment, or variance in expression, or effect size estimates from the model (Deviance), depending on what you want
                                                                                                   # sometimes here I like to save a figure conditional on there being a significant effect in one of the above
                                                                                                   if( <test whether Pop, Treat, or PopTreat P vales < 0.05/#_genes >){
                                                                                                       #add code to make a figure of your choice showing expression ~ population or treatment or both. Save the figure as a png named by the gene name
                                                                                                       # This is handy to go and visually check effect sizes, as a reality check that the P value you record actually “looks” plausible and isn’t caused by some weird outlier with large leverage
                                            }
                                            # end loop
                                            }
# name the results.table columns
names(results.table) <- c(“gene_name”, “mean_expression_relative”,  “Population_P”, “Treatment_P”, “PopTreatment_P"

# sort results table to highlight significant cases at the top
results.table <- results.table[order(results.table$PopTreatment_P, decreasing = F) , ]

# Show top 10 genes with most significant
results.table[1:10, ]