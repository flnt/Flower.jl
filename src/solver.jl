function init_ksp_solver(A, n, nullspace=false, ns_vec=nothing;
                         ksp_monitor_true_residual = false,
                         ksp_converged_reason = false,
                         ksp_constant_null_space = false)
    mat = PETSc.MatSeqAIJ(A)
    
    if nullspace && !isnothing(ns_vec)
        vecseq = PETSc.VecSeq(ns_vec)
        comm = PETSc.getcomm(mat)

        ns = PETSc.MatNullSpace{Float64}(comm, PETSc.PETSC_FALSE, 1, [vecseq])
        PETSc.MatSetNullSpace!(mat, ns)

        if n <= 600
            ksp = PETSc.KSP(mat;
                            ksp_monitor_true_residual = ksp_monitor_true_residual,
                            ksp_converged_reason = ksp_converged_reason,
                            ksp_constant_null_space = ksp_constant_null_space,
                            # ksp_type = "preonly",
                            # ksp_rtol = 1e-8,
                            # pc_type = "lu",
                            ksp_type = "cg",
                            ksp_rtol = 1e-8,
                            pc_type = "gamg",
                            pc_gamg_agg_nsmooths = 1,
                            pc_gamg_threshold = 0.02,
                            mg_levels_ksp_type = "chebyshev",
                            mg_levels_ksp_max_it = 5,
                            mg_levels_pc_type = "sor",
            )
        else
            ksp = PETSc.KSP(mat;
                            ksp_monitor_true_residual = ksp_monitor_true_residual,
                            ksp_converged_reason = ksp_converged_reason,
                            ksp_constant_null_space = ksp_constant_null_space,
                            ksp_type = "cg",
                            ksp_rtol = 1e-8,
                            pc_type = "gamg",
                            pc_gamg_agg_nsmooths = 1,
                            pc_gamg_threshold = 0.02,
                            mg_levels_ksp_type = "chebyshev",
                            mg_levels_ksp_max_it = 5,
                            mg_levels_pc_type = "sor",
            )
        end

        PETSc.destroy(vecseq)
        PETSc.destroy(mat)

        return ksp, ns
    else
        if n <= 600
            ksp = PETSc.KSP(mat;
                            ksp_monitor_true_residual = ksp_monitor_true_residual,
                            ksp_converged_reason = ksp_converged_reason,
                            ksp_constant_null_space = ksp_constant_null_space,
                            ksp_type = "preonly",
                            ksp_rtol = 1e-8,
                            pc_type = "lu",
            )
        else
            ksp = PETSc.KSP(mat;
                            ksp_monitor_true_residual = ksp_monitor_true_residual,
                            ksp_converged_reason = ksp_converged_reason,
                            ksp_constant_null_space = ksp_constant_null_space,
                            ksp_type = "cg",
                            ksp_rtol = 1e-8,
                            pc_type = "gamg",
                            pc_gamg_agg_nsmooths = 1,
                            pc_gamg_threshold = 0.02,
                            mg_levels_ksp_type = "chebyshev",
                            mg_levels_ksp_max_it = 5,
                            mg_levels_pc_type = "sor",
            )
        end

        PETSc.destroy(mat)

        return ksp, false
    end
end

function update_ksp_solver!(ksp, A, nullspace=false, ns_vec=nothing)
    mat = PETSc.MatSeqAIJ(A)

    if nullspace && !isnothing(ns_vec)
        vecseq = PETSc.VecSeq(ns_vec)
        comm = PETSc.getcomm(mat)

        ns = PETSc.MatNullSpace{Float64}(comm, PETSc.PETSC_FALSE, 1, [vecseq])
        PETSc.MatSetNullSpace!(mat, ns)

        PETSc.setoperators!(ksp, mat)
        PETSc.setfromoptions!(ksp)

        PETSc.destroy(vecseq)
        PETSc.destroy(mat)
        
        return ns
    else
        PETSc.setoperators!(ksp, mat)
        PETSc.setfromoptions!(ksp)

        PETSc.destroy(mat)

        return false
    end
end