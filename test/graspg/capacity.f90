! Check sizes across separately compiled suite and CSFG modules.
program csfg_capacity
    use parameter_def, only: NNNP, NNN1, NNNW
    use grid_C, only: R
    use def_C, only: DP
    use csfg_parameters
    use csfg_scf_C, only: UCF, METHOD, JXIPOS
    use csfg_tv_C, only: NBPT, NBCRPT, NBPV, NBCRPV, NCOLBLK, MCPDF1
    use symexpand_mod, only: TotCSFs_perblock, scfDF1, strlaborb
    use symmatrix_mod, only: NCFBLKTr, matrix_cutoff => CUTOFF
    use symmatrix_restart_C, only: restart_cutoff => CUTOFF
    use csfg_decide_C, only: CRITE0, NDISBR, LDISBR
    implicit none
    if (size(R) /= NNN1 .or. size(DP) /= NNNP) error stop 'radial capacity mismatch'
    if (size(UCF) /= NNNW .or. size(METHOD) /= NNNW) error stop 'orbital capacity mismatch'
    if (any(shape(JXIPOS) /= [NNNW+1, NNNW])) error stop 'orbital index shape mismatch'
    if (size(NCFBLKTr) /= CSFG_MAX_BLOCKS) error stop 'expanded block count capacity mismatch'
    if (any(shape(TotCSFs_perblock) /= [3, CSFG_MAX_BLOCKS])) error stop 'block table mismatch'
    if (size(scfDF1) /= CSFG_MAX_BLOCKS .or. size(MCPDF1) /= CSFG_MAX_BLOCKS) error stop 'DF1 size mismatch'
    if (any([size(NBPT), size(NBCRPT), size(NBPV), size(NBCRPV), size(NCOLBLK)] /= &
        CSFG_MAX_BLOCKS+1)) error stop 'block offset size mismatch'
    if (len(strlaborb) /= CSFG_LABEL_LENGTH) error stop 'label capacity mismatch'
    if (matrix_cutoff /= CSFG_MATRIX_CUTOFF .or. restart_cutoff /= CSFG_RESTART_CUTOFF) &
        error stop 'matrix cutoff mismatch'
    if (CRITE0 /= CSFG_DEFAULT_CRITE0 .or. NDISBR /= CSFG_DEFAULT_NDISBR .or. &
        LDISBR /= CSFG_DEFAULT_LDISBR) error stop 'default mismatch'
    print *, 'Shared CSFG capacities:', NNNP, NNN1, NNNW, CSFG_MAX_BLOCKS
end program csfg_capacity
