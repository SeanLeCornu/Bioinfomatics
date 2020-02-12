-- sean_exports.sql
-- Oct '19
-- Oracle based SQL for data exports
-- NB includes some "Old" (dt11c) data via legacy tables in new Ora db

    -- exp01[see exp04]: exp01_tcga_luad_lusc_top_mut_genes_ffs_20191003
        -- ** incomplete varmap data at this run - see expo04**
        SELECT
            mcf.data_source,
            tcga.cancer_type,
            mcf.sf_id,
            mcf.ff_id,
            mcf.rep_id,
            mcf.rep_source_id,
            mcf.vm_gene,            
            mcf.vm_uniprot_accession,
            mcf.mm_uniprot,
            mcf.variant_type,
            mcf.variant_class,
            mcf.vm_synonymous,
            count(*) mut_count_by_cancer_ff_gene
            
        FROM
            funvar_tx.mvw_map_core_funfam_pdb mcf
        LEFT JOIN FUNVAR_IMPORT.mut_tcga tcga
            ON mcf.hid = tcga.tcga_hid 
        where cancer_type in ('LUAD', 'LUSC')
            and mcf.variant_class= 'Missense_Mutation'
        group by mcf.data_source, tcga.cancer_type, mcf.sf_id, mcf.ff_id, mcf.rep_id, 
            mcf.rep_source_id, mcf.vm_gene, mcf.vm_uniprot_accession, mcf.mm_uniprot, mcf.variant_type, 
            mcf.variant_class, mcf.vm_synonymous
        having count(*) > 9
        ORDER BY tcga.cancer_type, sf_id,ff_id, rep_id, COUNT(*) DESC;
    
    -- exp02: PSEs (= FVEs "FunVar Events" now!) with full info for LUAD and LUSC
        -- old data / dt11c
        SELECT
            cancer_type,
            timing,
            LISTAGG(kegg_metab, '; ') WITHIN GROUP (ORDER BY kegg_metab) kegg_metabolic_pathways,
            gene,
            uniprot_accession,
            mut_pdb_res,
            mut_aa_chg,
            pse_count,
            tot_lung_muts,
            tot_pancan_muts,
            tot_pancan_muts_by_aachg,
            csa,
            ibis_small,
            ibis_nucl,
            ibis_ppi,
            scorecons_90,
            ff_dops,
            superfamily_id,
            funfam_number,
            rep_id,
            ff_name,
            neodriver_score_max_g_by_ff,
            tx_mutated_cgc_genes,
            ec3_num,
            ec4_num
        FROM
            funvar_tx.dt11c_dte_014_counts_pse 
        
        -- (exp02b)
            -- 21/10/19 - note can test what PSEs present in specific pathways enriched from Metascape, e.g. for https://reactome.org/PathwayBrowser/#/R-HSA-5579029&DTAB=MT
            --WHERE uniprot_accession in ('P22309','Q16678','P48637','Q03154','Q6V0L0','Q6NT55','P51580','O14841','P48507','P48506','Q02318','P23526','P31513','P22310','P15538','P19440','Q9NTN3','P08686','Q07973','Q7Z449','P19099','O15528','Q9NR63','P21397','O75881','P24557','P10109','Q6P4F2','P05108','P22570','Q6VVX0','Q00266','P11511','P01189','Q01718','P05093')

        group by 
            cancer_type, timing, gene, uniprot_accession, mut_pdb_res, 
            mut_aa_chg, pse_count, tot_lung_muts,tot_pancan_muts,
            tot_pancan_muts_by_aachg, csa, ibis_small, 
            ibis_nucl, ibis_ppi, scorecons_90, ff_dops, superfamily_id, funfam_number, 
            rep_id, ff_name, neodriver_score_max_g_by_ff, tx_mutated_cgc_genes, ec3_num, 
            ec4_num

        order by gene;

    -- exp03: PSEs (= FVEs "FunVar Events" now!) with basic info for LUAD and LUSC
        -- old data / dt11c
        SELECT
            DISTINCT
            cancer_type,
            timing,
            gene,
            uniprot_accession,
            superfamily_id,
            funfam_number,
            rep_id,
            
            ff_name,
            ff_dops
        FROM
            funvar_tx.dt11c_dte_014_counts_pse 
        order by gene;

    -- exp04: Top mutated genes (10+ mutations) in LUAD and LUSC with FunFam info
        SELECT
            mcf.data_source,
            tcga.cancer_type,
            
            mcf.sf_id,
            mcf.ff_id,
            mcf.rep_id,
            mcf.rep_source_id,
            --mcf.hugo_symbol,
            mcf.vm_gene,
            mcf.vm_uniprot_accession,
            --mcf.vm_seq_no,--mcf.vm_uniprot_aa,--mcf.vm_aa_change,--mcf.vm_change_type,
            mcf.mm_uniprot,
            mcf.variant_type,
            mcf.variant_class,
            mcf.vm_synonymous,
            count(*) mut_count_by_cancer_ff_gene
            --mcf.mm_aa_chg,--mcf.mm_funfam_id,--mcf.mm_input_seq_id,--mcf.mm_input_seq_pos,--mcf.mm_input_res_aa,--mcf.mm_aln_pos,--mcf.mm_output_seq_id,--mcf.mm_output_seq_pos,--mcf.mm_output_seq_aa,--mcf.mm_output_pdb_res,--mcf.uniprot_ffrep_res_match,--mcf.mm_info_text,--mcf.hid,
        --    mcf.mutation_id,--    mcf.data_import_id,--    mcf.map_import_id,--    mcf.mm_import_id,--    mcf.vm_note,--    mcf.ffr_member_id,--    mcf.ffr_aa_ranges,--    mcf.ffr_gene_name,--    mcf.ffr_import_id,--    mcf.mm_mutation_id,--    mcf.mm_hid
        FROM
            funvar_tx.mvw_map_core_funfam_pdb mcf

        LEFT JOIN FUNVAR_IMPORT.mut_tcga tcga
            ON mcf.hid = tcga.tcga_hid 
            
        where cancer_type in ('LUAD', 'LUSC')
        AND data_source = 'TCGA'
        and mcf.variant_class= 'Missense_Mutation'
            group by mcf.data_source, tcga.cancer_type, mcf.sf_id, mcf.ff_id, mcf.rep_id, 
        mcf.rep_source_id, mcf.vm_gene, mcf.vm_uniprot_accession, mcf.mm_uniprot, mcf.variant_type, 
        mcf.variant_class, mcf.vm_synonymous

        having count(*) > 9
        order by tcga.cancer_type, rep_id, count(*) desc;
    
    -- exp05: exp05_tcga_luad_lusc_top_mut_genes_ALL_FF_20191011
        -- ALL funfams (i.e. struct or seq)
        SELECT
            mcf.data_source,
            tcga.cancer_type,
            mcf.superfamily_id,
            mcf.funfam_number,
            mcf.rep_id,
            mcf.rep_source_id,
            mcf.vm_gene,            
            mcf.vm_uniprot_accession,
            mcf.variant_type,
            mcf.variant_class,
            mcf.vm_synonymous,
            count(*) mut_count_by_cancer_ff_gene
            
        FROM
            funvar_tx.mvw_map_core_funfam mcf
        LEFT JOIN FUNVAR_IMPORT.mut_tcga tcga
            ON mcf.hid = tcga.tcga_hid 
        where cancer_type in ('LUAD', 'LUSC')
            and mcf.variant_class= 'Missense_Mutation' 
            
        group by mcf.data_source, tcga.cancer_type, mcf.superfamily_id, mcf.funfam_number, mcf.rep_id, 
        mcf.rep_source_id, mcf.vm_gene, mcf.vm_uniprot_accession, mcf.variant_type, mcf.variant_class, 
        mcf.vm_synonymous
       
        having count(*) > 9
        ORDER BY tcga.cancer_type, superfamily_id,funfam_number, rep_id, COUNT(*) DESC;       
        
    -- exp06: exp06_tcga_luad_lusc_top_mut_genes_20191011.xlsx
        -- Top mut TCGA genes, without FunFam mapping
        -- Note count(*) threshold of 10+ here is effectively less stringent than the
        -- CATH FunFam queries (exp05, 04 etc) so could be raised when doing GSEAs
        SELECT
            mcf.data_source,
            tcga.cancer_type,
            mcf.vm_gene,            
            mcf.vm_uniprot_accession,
            mcf.variant_type,
            mcf.variant_class,
            mcf.vm_synonymous,
            count(*) mut_count_by_cancer_ff_gene
            
        FROM
            funvar_tx.mvw_map_core mcf
        LEFT JOIN FUNVAR_IMPORT.mut_tcga tcga
            ON mcf.hid = tcga.tcga_hid 
        where cancer_type in ('LUAD', 'LUSC')
            and mcf.variant_class= 'Missense_Mutation'
        
        group by 
            mcf.data_source, tcga.cancer_type, mcf.vm_gene, mcf.vm_uniprot_accession, mcf.variant_type, 
            mcf.variant_class, mcf.vm_synonymous 

        having count(*) > 9
        ORDER BY tcga.cancer_type, COUNT(*) DESC;    
        
    --- exp07 Background gene sets
        -- (a) All human FunFam genes seq or struct
        -- Gene name seems more accurate using UniProt core ( gene_names_primary )
        -- esp where there are multiple gene names per uniprot
        -- Check with HUGO if nec. e.g. 
        -- https://www.genenames.org/tools/search/#!/all?query=P30049
        --P30049	ATP5D	ATPD_HUMAN	ATP5F1D
        SELECT
                DISTINCT
                fug.uniprot_acc,
                fug.gene_name,
                uc.entry_name,
                uc.gene_names_primary
            FROM
                funvar_import.cath_funfam_to_uniprot_gene fug
            INNER JOIN FUNVAR_IMPORT.uniprot_core uc
                ON fug.uniprot_acc = uc.entry
            WHERE fug.taxon_id = 9606 and uc.status='reviewed'
                AND uc.gene_names_primary <> '9606'
                --AND fug.gene_name <> uc.gene_names_primary
            ORDER BY GENE_NAME;

    --- exp08: NFE core (display)
        -- 6th Nov 2019
        -- Generates core NFE fields as per:
        --    exp08_NFE_snapshot_20191106_for_imported_missense_mutclusts_plus_syn_part_fsite_20191106.xlsx
        --  **** Data will change as more MutClust FunFams imported and Functional site types added ****
        SELECT
            group_id NFE_group,
            groupid mutclut_group,
            sf_id,
            ff_id,
            rep_id,
            rep_source_id,
            hugos_t,
            --uniprots_t,
            --seq_nos_t,
            aa_changes_t,
            gonmad_afs_t,
            nat_variants_t,

            cadd_marks_t,
            mutclust_residue,
            on_scons_90,
            on_ppi,
            on_nuc,
            on_lig,
            on_lig_id,
            near_nuc,
            near_nuc_dist,
            near_lig,
            near_lig_dist,
            near_lig_id,
            diseases_t,
            diseae_variants_t,
            groupid,
            runid,
            taskid,
            funfam_name,
            mfc_mut_count,
            num_swissprot_in_ff,
            num_cgc_genes,
            atom_length,
            row_count,
            sig_mut_count_sum_p_corr,
            sig_weighted_mut_sum_p_corr,
            mut_count_sum_p,
            mut_count_sum_p_corr,
            mut_count_sum_p_corr_sig,
            weighted_mut_sum_p,
            weighted_mut_sum_p_corr,
            weighted_mut_sum_p_corr_sig,
            sc90_fsite_type,
            sc90_radius,
            sc90_ref_pdb_label,
            sc90_ref_pdb_near,
            ppi_fsite_type,
            ppi_radius,
            ppi_ref_pdb_label,
            ppi_ref_pdb_near,
            nuc_fsite_type,
            nuc_radius,
            nuc_ref_pdb_label,
            nuc_ref_pdb_near,
            lig_fsite_type,
            lig_radius,
            lig_ligand,
            lig_ref_pdb_label,
            lig_ref_pdb_near,
            near_nuc_fsite_type,
            near_nuc_radius,
            near_nuc_ref_pdb_label,
            near_nuc_ref_pdb_near,
            near_lig_fsite_type,
            near_lig_radius,
            near_ligand,
            near_lig_ref_pdb_label,
            near_lig_ref_pdb_near,
            mm_output_pdb_res,
            mutids
        FROM
            funvar_tx.vw_nfe_core_display
        WHERE
            -- Only show if have missense MutClusts returned & imported
            taskid IN(
                SELECT DISTINCT taskid FROM funvar_import.mutclust_analysis WHERE groupid='mc07')

        ORDER BY mfc_mut_count DESC, weighted_mut_sum_p_corr;

    -- exp09: NFE v3 / 2 
        -- 21/11/2019
        -- NFE mutations with mutclust runs to ~19/11
        -- Includes cancer type, GD info and data source
        -- Only ON func site NFEs (no near site data in base tables)
        -- Be aware that one mutation_id (NFEs here) can have many ligand_ids, if changing SELECT list.
        SELECT
            nfe_version,
            data_source,
            cancer_type,
            sf_id,
            ff_id,
            rep_id,
            funfam_name,
            hugo_symbol,
            
            --  NOT distinct on mutation_ids because of ligand IDs - use DISTINCT in related queries e.g.
            COUNT( DISTINCT mutation_id ) num_NFE_per_ff_gene,            
            
            vm_uniprot_accession,
            mc_num_muts,
            mc_num_swissprot_ff,
            mc_num_cgc_genes,
            --mc_num_res_ff_rep,
            --mc_num_mut_res_ff_rep,
            mc_sig_mut_count_sum_corr,
            mc_sig_weighted_mut_sum_corr
        
        FROM
            funvar_archive.nfe_03

        WHERE 
            cancer_type IN ('LUAD', 'LUSC')
            --cancer_type ='LUAD'

            -- MISSENSE, on or near any f-site
            -- Current live MutClust for missense muts
            AND groupid = 'mc07'
                -- Make sure on/near any functional site
                AND NOT
                (
                    -- Mutclust mutation is not ON any functional site
                    -- NB -1: is for sites not yet implemented
                        ( on_scons_90 = 0 OR on_scons_90 = -1 )
                    AND ( on_mcsa = 0 OR on_mcsa = -1 )  
                    AND ( on_ppi = 0 OR on_ppi = -1 )
                    AND ( on_nuc = 0 OR on_nuc = -1 )
                    AND ( on_lig = 0 OR on_lig = -1 )

                    -- Mutclust mutation is not NEAR any functional site
                    -- NB -1: is for sites not yet implemented
                    AND ( near_scons_90 = 0 OR near_scons_90 = -1 )
                    AND ( near_mcsa = 0  OR near_mcsa = -1 )
                    AND ( near_nuc = 0 OR near_nuc = -1 )
                    AND ( near_lig = 0 OR near_lig = -1 )
                
                ) 
                
            -- Optional 
            -- AND gdstatus = 'GD'   -- GD tumours only
            -- AND gdstatus = 'nGD'   -- non GD tumours only
            -- AND timing = 'early'
            -- etc
            
        GROUP BY 
            nfe_version, data_source, cancer_type, sf_id, ff_id, 
            rep_id, funfam_name, hugo_symbol, vm_uniprot_accession, mc_num_muts, 
            mc_num_swissprot_ff, mc_num_cgc_genes, mc_sig_mut_count_sum_corr, mc_sig_weighted_mut_sum_corr
        
        ORDER BY cancer_type, hugo_symbol;
        
    
    -- exp10 - copy of SQL used in neofun/script/sql/analysis/NFE_analysis.sql to give summary counts 
    -- in "NFE_Figures.key" for nfe_v04
    -- Comment WHERE clause in/out as appropriate!
    -- 09/12/19
     -------------------
	-- NFE snapshot v04
	-------------------
    -- FOR TABLE
        -- Mutations (MS)
            SELECT COUNT(*) FROM
                ( SELECT distinct mutation_id FROM funvar_archive.nfe_04
                    WHERE groupid='mc07' 
                    --AND cancer_type='LUAD'
                    --AND cancer_type='LUSC'
                    );
                
        -- NFE-ON
            SELECT COUNT(*) FROM
                ( SELECT distinct mutation_id FROM funvar_archive.nfe_04
                    WHERE groupid='mc07' 
                        --AND cancer_type='LUAD'
                    --AND cancer_type='LUSC'
                        AND NOT
                        (   ( on_scons_90 = 0 OR on_scons_90 = -1 ) AND ( on_mcsa = 0 OR on_mcsa = -1 )  AND 
                            ( on_ppi = 0 OR on_ppi = -1 ) AND ( on_nuc = 0 OR on_nuc = -1 ) AND ( on_lig = 0 OR on_lig = -1 )
                        ) 
                );

		-- Genes
                SELECT COUNT(*) FROM
                    ( SELECT distinct hugo_symbol FROM funvar_archive.nfe_04
                        WHERE groupid='mc07' 
                        --AND cancer_type='LUAD'
                       -- AND cancer_type='LUSC'
                        );
                    
			-- NFE-ON
                SELECT COUNT(*) FROM
                    ( SELECT distinct hugo_symbol FROM funvar_archive.nfe_04
                        WHERE groupid='mc07' 
                         --AND cancer_type='LUAD'
                        --AND cancer_type='LUSC'
                            AND NOT
                            (   ( on_scons_90 = 0 OR on_scons_90 = -1 ) AND ( on_mcsa = 0 OR on_mcsa = -1 )  AND 
                                ( on_ppi = 0 OR on_ppi = -1 ) AND ( on_nuc = 0 OR on_nuc = -1 ) AND ( on_lig = 0 OR on_lig = -1 )
                            ) 
                    );

		-- FunFams
                SELECT COUNT(*) FROM
                    ( SELECT distinct rep_id FROM funvar_archive.nfe_04
                        WHERE groupid='mc07' 
                        --AND cancer_type='LUAD'
                        --AND cancer_type='LUSC'
                        );
                    
			-- FVE-ON
                SELECT COUNT(*) FROM
                    ( SELECT distinct rep_id FROM funvar_archive.nfe_04
                        WHERE groupid='mc07' 
                         --AND cancer_type='LUAD'
                        --AND cancer_type='LUSC'
                            AND NOT
                            (   ( on_scons_90 = 0 OR on_scons_90 = -1 ) AND ( on_mcsa = 0 OR on_mcsa = -1 )  AND 
                                ( on_ppi = 0 OR on_ppi = -1 ) AND ( on_nuc = 0 OR on_nuc = -1 ) AND ( on_lig = 0 OR on_lig = -1 )
                            ) 
                    );
		
		-- NFEsites (i.e. single positions / mutclusts on FunFam rep)
            -- TOTAL
            -- TOTAL_W 
                SELECT COUNT(*) FROM
                    ( SELECT distinct rep_id, mutclust_residue FROM funvar_archive.nfe_04
                        WHERE groupid='mc07' 
                        --AND cancer_type='LUAD'
                        --AND cancer_type='LUSC'
                        --AND weighted_mut_sum_p_corr_sig ='Y'
                        );
                    
			-- NFE-ON
            -- NFE-ON_W 
                SELECT COUNT(*) FROM
                    ( SELECT distinct rep_id, mutclust_residue FROM funvar_archive.nfe_04
                        WHERE groupid='mc07' 
                        AND cancer_type='LUAD'
                        --AND cancer_type='LUSC'
                        -- AND weighted_mut_sum_p_corr_sig ='Y'
                            AND NOT
                            (   ( on_scons_90 = 0 OR on_scons_90 = -1 ) AND ( on_mcsa = 0 OR on_mcsa = -1 )  AND 
                                ( on_ppi = 0 OR on_ppi = -1 ) AND ( on_nuc = 0 OR on_nuc = -1 ) AND ( on_lig = 0 OR on_lig = -1 )
                            ) 
                    );


        -- exp11 : NFE_v06 "beta" (using live materialized view - this could change, unlike data in funvar_archive)
        -- Lung cancer NFEs (NFE + PFH) for mis-sense muts
        -- NOT FunVar scored; PFHs have same amino acid change, same gene in more than 3+ patients
        -- NOTE 23/01/20: erroneously commented out "on_lig" and "near_lig" 
        --  these can be included for testing if NFE ON or NEAR ligand site!
            SELECT DISTINCT
                mutation_id, hid,
                nfe_type, nfe_version, data_source, cancer_type, variant_type, variant_class, vm_synonymous, 
                sf_id, ff_id, rep_id, mutfam, mutfam_version, mfc_mut_count_missense, mfc_mut_count_silent, 
                hugo_symbol, vm_uniprot_accession, vm_seq_no, vm_aa_change, pdb_res, 
                on_scons_90, on_mcsa, on_ppi, on_nuc, -- on_lig, on_lig_id, 
                near_angstroms, near_scons_90, near_mcsa, near_nuc, -- near_lig, near_lig_id, 
                num_gd, num_ngd, num_timing_early, num_timing_late, num_timing_not_poss, clonality, num_clonal, num_subclonal, 
                gnomad_af, diseases, cadd_mark, protein_cat_res, protein_disulphide, mc_num_muts, num_patients, 
                mut_count_sum_p, mut_count_sum_p_corr, mut_count_sum_p_corr_sig, 
                weighted_mut_sum_p, weighted_mut_sum_p_corr, weighted_mut_sum_p_corr_sig, groupid, runid, taskid, 
                mc_num_swissprot_ff, mc_num_cgc_genes, mc_num_res_ff_rep, mc_num_mut_res_ff_rep, mc_sig_mut_count_sum_corr, mc_sig_weighted_mut_sum_corr
                FROM
                    funvar_tx.mvw_nfe_pfh_mutfam
                WHERE
                    vm_synonymous='FALSE' 
                    AND ( ( nfe_type = 'MC' AND groupid = 'mc07' ) OR ( nfe_type='PFH' AND num_patients > 2 ) )
                    AND cancer_type IN ('LUAD', 'LUSC') ;

        -- exp12 : Top ranked NFE score for each mutation position / AA change
        -- 28/01/2020
        -- This  has been made available as a VIEW; the view's SQL is here for reference!
        -- *** 11/02/2020 Note use the archive version as a fixed v6_1 ***
            ----SELECT * FROM funvar_tx.vw_nfe_pfh_scored_rank_mut;
            SELECT * FROM funvar_archive.nfe_pfh_scored_rank_mut_061;
            
            --DROP VIEW funvar_tx.vw_nfe_pfh_scored_rank_mut;
            -- CREATE VIEW funvar_tx.vw_nfe_pfh_scored_rank_mut
            -- AS
            -- WITH score_ranked AS
            -- (
            -- SELECT
            --     nfe_version,
            --     data_source,
            --     cancer_type,
            --     nfe_type,
            --     variant_type,
            --     variant_class,
            --     vm_synonymous,
            --     sf_id,
            --     ff_id,
            --     rep_id,
            --     mutfam,
            --     hugo_symbol,
            --     vm_uniprot_accession,
            --     vm_seq_no,
            --     pdb_res,
            --     vm_aa_change,
            --     hotspot_num_muts,
            --     nfe_score_d_mf_h,
            --     gnomad_af,
            --     ispoly,
            --     deltasize,
            --     mclachlan,
            --     diseases,

            --     on_scons_90,
            --     on_mcsa,
            --     on_ppi,
            --     on_nuc,
            --     on_lig,
            --     near_angstroms,
            --     near_scons_90,
            --     near_mcsa,
            --     near_nuc,
            --     near_lig,
            --     cadd_mark,
            --     num_gd,
            --     num_ngd,
            --     num_timing_early,
            --     num_timing_late,
            --     num_timing_not_poss,   
                
            --     ROW_NUMBER() OVER (
            --             PARTITION BY cancer_type, sf_id, ff_id, rep_id, hugo_symbol, vm_seq_no, vm_aa_change 
            --             ORDER BY cancer_type, sf_id, ff_id, rep_id, hugo_symbol, vm_seq_no, vm_aa_change, nfe_score_d_mf_h DESC 
            --         ) AS rn
    
            --     FROM
            --         funvar_tx.vw_nfe_pfh_scored
            --         ORDER BY cancer_type, sf_id, ff_id, rep_id, hugo_symbol, vm_seq_no, vm_aa_change
            --     )
            --     SELECT *
            --     FROM
            --         score_ranked 
            --     WHERE
            --         rn = 1 ;
    
	-- exp13 : Retrieve TCGA barcodes
		SELECT DISTINCT
    			TCGABARCODE from FUNVAR_IMPORT.mut_tcga where cancer_type = 'LUAD' ;


    -- exp14: Scored v6_1 NFEs with focus on the GD and non-GD fields
        SELECT
            nfe_version,
            data_source,
            cancer_type,
            --variant_type,
            --variant_class,
            --vm_synonymous,
            sf_id, ff_id, rep_id, mutfam, 
            hugo_symbol, vm_uniprot_accession, vm_seq_no, pdb_res, vm_aa_change, nfe_score_d_mf_h, 
            nfe_type, hotspot_num_muts, num_gd, num_ngd, num_timing_early, num_timing_late, num_timing_not_poss,
            on_scons_90, on_mcsa, on_ppi, on_nuc, on_lig, 
            near_angstroms, near_scons_90, near_mcsa, near_nuc, near_lig, 
            cadd_mark, gnomad_af, ispoly, deltasize, mclachlan, diseases, rn, mutation_id
        FROM
            funvar_archive.nfe_pfh_scored_rank_mut_061
        --WHERE cancer_type='LUSC'
        ORDER BY cancer_type,rep_id, hugo_symbol, nfe_score_d_mf_h DESC;
