function sparams = solveCMEsSingleParticlefromSameOrbital( sparams )
%SOLVECMESSINGLEPARTICLEFROMSAMEORBITAL Summary of this function goes here
%   Detailed explanation goes here
    % In an effort to reduce memory overhead, we will not include spin for
    % now (since there is no magnetic field this is easy to add in later).
    sparams.CMEsSingleParticle = zeros(sparams.nNonShiftedHOs^2, sparams.nNonShiftedHOs^2);
    
    as = sparams.acoeffs;
    
    % This part is fairly simple in code.
    % We wish to convert our CMEs that we found in the non-shifted HO basis
    % to the basis of the single particle orbitals for our arbitrary
    % potentials.  We have two matrices to help us out here:
    % A = sparams.acoeffs and B = sparams.bcoeffs.  A transforms the
    % single particle orbitals to the shifted HO states.  B tranforms the
    % shifted HO states to the non-shifted HO basis states.
    
    % First, we will transform the non-shifted CMEs to the shifted HO
    % basis.
    % To understand why we do what we do here, consider that the matrix of
    % CMEs in the non-shifted basis looks like:
    % [<a_1b_1|v|g_1d_1>, <a_1b_1|v|g_1d_2> ... <a_1b_1|v|g_1d_M>, <a_1b_1|v|g_2d_1> ... <a_1b_1|v|g_Md_M>]
    % [<a_1b_2|v|g_1d_1>, <a_1b_2|v|g_1d_2> ... <a_1b_2|v|g_1d_M>, <a_1b_2|v|g_2d_1> ... <a_1b_2|v|g_Md_M>]
    % [       ...                ...                   ...                ...                   ...       ]
    % [<a_1b_M|v|g_1d_1>, <a_1b_M|v|g_1d_2> ... <a_1b_M|v|g_1d_M>, <a_1b_M|v|g_2d_1> ... <a_1b_M|v|g_Md_M>]
    % [<a_2b_1|v|g_1d_1>, <a_2b_1|v|g_1d_2> ... <a_2b_1|v|g_1d_M>, <a_2b_1|v|g_2d_1> ... <a_2b_1|v|g_Md_M>]
    % [       ...                ...                   ...                ...                   ...       ]
    % [<a_Mb_M|v|g_1d_1>, <a_Mb_M|v|g_1d_2> ... <a_Mb_M|v|g_1d_M>, <a_Mb_M|v|g_2d_1> ... <a_Mb_M|v|g_Md_M>]
    
    % The B matrix looks like:
    % [<r_1|a_1>, <r_1|a_2> ... <r_1|a_M>]
    % [<r_2|a_1>, <r_2|a_2> ... <r_2|a_M>]
    % [   ...        ...           ...   ]
    % [<r_N|a_1>, <r_N|a_2> ... <r_N|a_M>]
    
    % When we tensor B x B we obtain the matrix:
    % [<r_1s_1|v|a_1b_1>, <r_1s_1|v|a_1b_2> ... <r_1s_1|v|a_1b_M>, <r_1s_1|v|a_2b_1> ... <r_1s_1|v|a_Mb_M>]
    % [<r_1s_2|v|a_1b_1>, <r_1s_2|v|a_1b_2> ... <r_1s_2|v|a_1b_M>, <r_1s_2|v|a_2b_1> ... <r_1s_2|v|a_Mb_M>]
    % [       ...                ...                   ...                ...                   ...       ]
    % [<r_1s_N|v|a_1b_1>, <r_1s_N|v|a_1b_2> ... <r_1s_N|v|a_1b_M>, <r_1s_N|v|a_2b_1> ... <r_1s_N|v|a_Mb_M>]
    % [<r_2s_1|v|a_1b_1>, <r_2s_1|v|a_1b_2> ... <r_2s_1|v|a_1b_M>, <r_2s_1|v|a_2b_1> ... <r_2s_1|v|a_Mb_M>]
    % [       ...                ...                   ...                ...                   ...       ]
    % [<r_Ns_N|v|a_1b_1>, <r_Ns_N|v|a_1b_2> ... <r_Ns_N|v|a_1b_M>, <r_Ns_N|v|a_2b_1> ... <r_Ns_N|v|a_Mb_M>]
    
    % Written out this way, it is clear to see that in order to convert the
    % CMEs to the {r,s,t,u} basis (localized HO basis), we simply do:
    % (B^(-1) x B^(-1))*CME*(B x B) where CME is the non-shifted CME matrix
    % Second, we will transform the shifted HO basis to the single particle
    % orbital basis using the sample procedure but with the acoeffs matrix
    % instead

%     sparams.CMEsSingleParticle = kron(inv(as),inv(as))*sparams.CMEsLOHO*kron(as,as);
    sparams.CMEsSingleParticle = kron(as',as')*sparams.CMEsNonShifted*kron(as,as);
end

