{
    labelList cellBoundaryFaceCount(omega_.size(), 0);
    
    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                omega_[faceCelli] = 0.0;
            }
        }
    }
    
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                cellBoundaryFaceCount[faceCelli]++;

                omega_[faceCelli] =
                    2.0*this->nu()()[faceCelli]/sqr(yr_[faceCelli]);
            }
        }
    }

    // perform the averaging

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                omega_[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }

    forAll(patches, patchi)
    {
        const fvPatch& p = patches[patchi];

        if (isType<wallFvPatch>(p))
        {
            omegaEqn->setValues
            (
                p.faceCells(),
                omega_.boundaryField()[patchi].patchInternalField()
            );
        }
    }
}
