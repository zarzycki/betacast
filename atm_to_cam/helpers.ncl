;************************************************************************* 

undef ("z_to_z_interp")
function z_to_z_interp (theta_fv:numeric,z_fv:numeric,thisCol:numeric,extrapLow:logical,extrapHigh:logical)
local nummiss, deriv, ff, nlev
begin

    nlev = dimsizes(thisCol)
    
    t_wrf = wrf_interp_1d (theta_fv,z_fv,thisCol)
    
    ; if we have any missing values, we need to extrapolate
    if any(ismissing(t_wrf)) then
      ;print(t_wrf(:,ix) + " " + thisCol)
      ;print("---")
      ixvalid = ind(.not.(ismissing(t_wrf)))

      if (ismissing(t_wrf(0))) then
        lowest_valid = min(ixvalid)
        
        ; if extrap is chosen, calculate true derivative, otherwise set to 0 for iso() atmosphere.
        if (extrapLow) then
          deriv = (t_wrf(lowest_valid) - t_wrf(lowest_valid+1)) / (thisCol(lowest_valid) - thisCol(lowest_valid+1))
        else
          deriv = 0.0
        end if
        
        do ff = 0,lowest_valid-1
          t_wrf(ff) = t_wrf(lowest_valid) + deriv * ( thisCol(ff) - thisCol(lowest_valid))
        end do
      end if
      
      if (ismissing(t_wrf(nlev-1))) then
        highest_valid = max(ixvalid)
        if (extrapHigh) then
          deriv = (t_wrf(highest_valid-1) - t_wrf(highest_valid)) / (thisCol(highest_valid-1) - thisCol(highest_valid))
        else
          deriv = 0.0
        end if
        do ff = highest_valid+1,nlev-1
          t_wrf(ff) = t_wrf(highest_valid) + deriv * ( thisCol(ff) - thisCol(highest_valid))
        end do
      end if
            
    end if
    
    return(t_wrf)
    
end



undef ("uv_cell_to_edge")
function uv_cell_to_edge (uZonal:numeric,uMerid:numeric,nlev:integer,lonEdge:numeric,latEdge:numeric,lonCell:numeric,latCell:numeric,edgeNormalVectors:numeric,cellsOnEdge:numeric)
local nCells
begin

    nCells = dimsizes(lonCell)
    nEdges = dimsizes(lonEdge)
    
    east = new((/3,nCells/),"double")
    east = 0.
    north = east
    
     do iCell = 0, nCells - 1

        east(0,iCell) = -sin(lonCell(iCell))
        east(1,iCell) =  cos(lonCell(iCell))
        east(2,iCell) =  0.0

        ; Normalize
        east(0:2,iCell) = east(0:2,iCell) / sqrt(sum(east(0:2,iCell) * east(0:2,iCell)))

        north(0,iCell) = -cos(lonCell(iCell))*sin(latCell(iCell))
        north(1,iCell) = -sin(lonCell(iCell))*sin(latCell(iCell))
        north(2,iCell) =  cos(latCell(iCell))

        ; Normalize
        north(0:2,iCell) = north(0:2,iCell) / sqrt(sum(north(0:2,iCell) * north(0:2,iCell)))

     end do
       
    uNormal = new((/nlev,nEdges/),typeof(uZonal))

     do iEdge = 0, nEdges - 1
        cell1 = cellsOnEdge(iEdge,0)-1
        cell2 = cellsOnEdge(iEdge,1)-1

        uNormal(:,iEdge) =  uZonal(:,cell1)      *0.5*(edgeNormalVectors(iEdge,0)*east(0,cell1)   \
                                                     + edgeNormalVectors(iEdge,1)*east(1,cell1)   \
                                                     + edgeNormalVectors(iEdge,2)*east(2,cell1))  \
                          + uMerid(:,cell1)      *0.5*(edgeNormalVectors(iEdge,0)*north(0,cell1)  \
                                                     + edgeNormalVectors(iEdge,1)*north(1,cell1)  \
                                                     + edgeNormalVectors(iEdge,2)*north(2,cell1)) \
                          + uZonal(:,cell2)      *0.5*(edgeNormalVectors(iEdge,0)*east(0,cell2)   \
                                                     + edgeNormalVectors(iEdge,1)*east(1,cell2)   \
                                                     + edgeNormalVectors(iEdge,2)*east(2,cell2))  \
                          + uMerid(:,cell2)      *0.5*(edgeNormalVectors(iEdge,0)*north(0,cell2)  \
                                                     + edgeNormalVectors(iEdge,1)*north(1,cell2)  \
                                                     + edgeNormalVectors(iEdge,2)*north(2,cell2))
     end do
     

return (uNormal)

end

