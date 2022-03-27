;************************************************************************* 

undef ("z_to_z_interp")
function z_to_z_interp (theta_fv:numeric,z_fv:numeric,thisCol:numeric,extrapLow:logical,extrapHigh:logical)
local nummiss, deriv, ff, nlev
begin

    nlev = dimsizes(thisCol)
    
    t_wrf = wrf_interp_1d (theta_fv,z_fv,thisCol)
    ;t_wrf = todouble(linint1(z_fv,theta_fv,False,thisCol,0))
    
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
        if (mod(iEdge,10000) .eq. 0) then
          print("... UV_CELL_TO_EDGE: "+(100.*iEdge/(nEdges-1))+"%")
        end if
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



;************************************************************************* 

undef ("add_time_define_precision")
function add_time_define_precision (var_in:numeric,precision:string,isncol:logical)
local ncol, nlev, mydims
begin

mydims = dimsizes(var_in)

if (isncol) then

  if (dimsizes(mydims) .eq. 1) then
  
    print("converting a ncol -> ncol "+precision+" with time attached!")
  
    ncol = mydims(0)
    if (precision .eq. "float" .or. precision .eq. "single") then
      var_out = new((/1,ncol/),float)
      var_out(0,:) = tofloat(var_in)
    else if (precision .eq. "double") then
      var_out = new((/1,ncol/),double)
      var_out(0,:) = todouble(var_in)
    else
      print("invalid precision, exiting...")
      exit
    end if
    end if
    
    var_out!0 = "time"
    var_out!1 = "ncol"
    
  else
  
    print("converting a lev,ncol -> lev,ncol "+precision+" with time attached!")

    nlev = mydims(0)
    ncol = mydims(1)
    
    if (precision .eq. "float" .or. precision .eq. "single") then
      var_out = new((/1,nlev,ncol/),float)
      var_out(0,:,:) = tofloat(var_in)
    else if (precision .eq. "double") then
      var_out = new((/1,nlev,ncol/),double)
      var_out(0,:,:) = todouble(var_in)
    else
      print("invalid precision, exiting...")
      exit
    end if
    end if
    
    var_out!0 = "time"
    var_out!1 = "lev"
    var_out!2 = "ncol"

  end if
  
  
else  ; if not ncol

  if (dimsizes(mydims) .eq. 2) then
  
    print("converting a nlat,nlon -> nlat,nlon "+precision+" with time attached!")
  
    nlat = mydims(0)
    nlon = mydims(1)
    if (precision .eq. "float" .or. precision .eq. "single") then
      var_out = new((/1,nlat,nlon/),float)
      var_out(0,:,:) = tofloat(var_in)
    else if (precision .eq. "double") then
      var_out = new((/1,nlat,nlon/),double)
      var_out(0,:,:) = todouble(var_in)
    else
      print("invalid precision, exiting...")
      exit
    end if
    end if
    
    var_out!0 = "time"
    var_out!1 = "lat"
    var_out!2 = "lon"
        
  else
  
    print("converting a lev,nlat,nlon -> lev,nlat,nlon "+precision+" with time attached!")

    nlev = mydims(0)
    nlat = mydims(1)
    nlon = mydims(2)
    
    if (precision .eq. "float" .or. precision .eq. "single") then
      var_out = new((/1,nlev,nlat,nlon/),float)
      var_out(0,:,:,:) = tofloat(var_in)
    else if (precision .eq. "double") then
      var_out = new((/1,nlev,nlat,nlon/),double)
      var_out(0,:,:,:) = todouble(var_in)
    else
      print("invalid precision, exiting...")
      exit
    end if
    end if
    
    var_out!0 = "time"
    var_out!1 = "lev"
    var_out!2 = "lat"
    var_out!3 = "lon"

  end if
  

end if

return(var_out)
    
end




;************************************************************************* 

undef ("latlon_to_ncol")
function latlon_to_ncol(var_in:numeric)
local nlev, nlat, nlon, ncol, vardims
begin

vardims = dimsizes(var_in)

if (dimsizes(vardims) .eq. 2) then
  var_out = ndtooned(var_in)
else if (dimsizes(vardims) .eq. 3) then
  nlev = vardims(0)
  nlat = vardims(1)
  nlon = vardims(2)
  ncol = nlat*nlon
  print("unpacking -> nlev: "+nlev+"    nlat: "+nlat+"    nlon: "+nlon+"    ncol: "+ncol)
  var_out = new((/nlev,ncol/),typeof(var_in))
  do ii = 0,nlev-1
    var_out(ii,:) = ndtooned(var_in(ii,:,:))
  end do
else
  print(vardims+" dims not supported")
  exit
end if
end if

return(var_out)
    
end


;************************************************************************* 

undef ("ncol_to_latlon")
function ncol_to_latlon(var_out:numeric,nlat:numeric,nlon:numeric)
local nlev, nlat, nlon, ncol, vardims
begin

vardims = dimsizes(var_out)

if (dimsizes(vardims) .eq. 1) then
  var_in = onedtond(var_out,(/nlat,nlon/))
else if (dimsizes(vardims) .eq. 2) then
  nlev = vardims(0)
  ncol = vardims(1)
  print("repacking -> nlev: "+nlev+"    nlat: "+nlat+"    nlon: "+nlon+"    ncol: "+ncol)
  var_in = new((/nlev,nlat,nlon/),typeof(var_out))
  do ii = 0,nlev-1
    var_in(ii,:,:) = onedtond(var_out(ii,:),(/nlat,nlon/))
  end do
else
  print(vardims+" dims not supported")
  exit
end if
end if

return(var_in)
    
end




