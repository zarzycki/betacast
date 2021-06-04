;************************************************************************* 
; D. Shea 
; Read simple f90 namelist file as indicated by Ben Foster [HAO]
; An integer variable is returned. Upon return it will have 
; associated with it attributes that were created from the
; input namelist file.
;
; Usage:
; 
;     hao_user = namelist("haoUserNameList")
;
;     In some cases, these new attributes are used to modify an
;     existing (say, default) suite of attributes (say, haoDefault).
;
;     hao = merge_VarAtts(hao_user, haoDefault) 
;
undef ("namelist")
function namelist (fname:string)
local lines, nam, chr, blankc, eqc, quot1c, quot2c, semic    \
    , decimalc, slashc, commac, ampc, nullc, newlc, nLines   \
    , nl, cc, CC, nc, NC, nEq, iEq, iSemic, aName, iq, i, n  \
    , iComma, nComma, nDecimal, rhs, RHS 
begin
  lines  = asciiread (fname, -1, "string")
  nam    = 1   ; create variable to which 
               ; attributes may be attached

               ; special characters to check 
               ; for clarity use stringtochar 
  chr    = stringtochar(" ")  ; dimsizes(chr)=2
  blankc = chr(0)             ; blank space     (int 32)

  chr    = stringtochar("=")  ; dimsizes(chr)=2
  eqc    = chr(0)             ; equal           (int 61)

  chr    = stringtochar(";")  ; dimsizes(chr)=2
  semic  = chr(0)             ; comment         (int 59)

  chr    = stringtochar("\")  ; dimsizes(chr)=2
  slashc = chr(0)             ; line continuation (int 92)

  chr    = stringtochar(",")  ; dimsizes(chr)=2
  commac = chr(0)             ; value separator (int 44)

  chr    = stringtochar(".")  ; dimsizes(chr)=2
  decimalc= chr(0)            ; indicate float  (int 46)

  chr    = stringtochar("&")  ; dimsizes(chr)=2
  ampc   = chr(0)             ; ampersand       (int 38)

  chr    = stringtochar("'")  ; dimsizes(chr)=2
  quot1c = chr(0)             ; single quote    (int 39)

  quot2c= integertochar(34)   ; double quote "  (int 34)
  newlc = integertochar(10)   ; new line character
  nullc = integertochar(0)    ; null

  nLines = dimsizes(lines)    ; # of lines (strings)

  do nl=0,nLines-1                 ; loop over each line
     cc  = stringtochar(lines(nl)) ; convert to characters
     nc  = dimsizes(cc)            ; # characters
     nc  = nc - 1                  ; ignore last character

     nEq = num(cc.eq.eqc)          ; number of = signs
                                   ; eliminate (a) HAO's  &proc
                                   ; (b) any line without an =
     if (cc(0).eq.ampc .or. nEq.eq.0) then     
         delete (cc)               ; delete cuz size changes
         continue                  ; go to next iteration
     end if

     iSemic = ind(cc.eq.semic)     ; is simicolon (;) present [comment] 
     if (.not.ismissing(iSemic)) then
         nc = iSemic(0)-1          ; only look up to semi-colon
         if (nc.le.1) then         ; need at least 3 char [eg: i=0]
             delete (cc)
             delete (iSemic)
             continue              ; go to next iteration
         end if
     end if
     delete (iSemic)

     NC = -1                       ; remove blanks
     CC = new ( nc, "character")   ; cc after blanks removed
   do n=0,nc-1
      if (cc(n).ne.blankc .and. cc(n).ne.nullc .and. cc(n).ne.newlc) then 
          NC     = NC+1
          CC(NC) = cc(n)
      end if
   end do
   delete (cc)                     ; no longer needed

   if (NC.le.1) then               ; again need min 3 char
       delete (CC)                 ; size might change
       continue                    ; go to next iteration
   end if

   iEq     = ind(CC.eq.eqc)        ; = is separator; return index

                                   ; name of attribute (lhs)
   aName    = chartostring( (/CC(0:iEq-1)/) )

   nComma   = num(CC.eq.commac)    ; a comma (,) means more than one 
   RHS      = chartostring(CC(iEq+1:))       ; right hand side

                              ; does rhs have a ' or " [if so, string]
   iq = ind(CC.eq.quot1c .or. CC.eq.quot2c)  ; indices of quotes (',")

   if (any(.not.ismissing(iq))) then
       CC(iq) = quot2c        ; change ' to " [also "-to-" for convenience]

       if (nComma.gt.0) then              ; more than 1 element
           rhs         = new ( nComma+1, "string")
           delete (rhs@_FillValue)
                                          ; must parse CC(iEq+1:)
                                          ; put each element => array rhs
           iComma      = ind(CC.eq.commac); indices of commas
           rhs(0)      = chartostring( CC(iEq+1:iComma(0)-1) )
           rhs(nComma) = chartostring( CC(iComma(nComma-1)+1:nc-1) )
           if (nComma.gt.1) then
               do i=0,nComma-2
                  rhs(i+1) = chartostring(CC(iComma(i)+1:iComma(i+1)-1))
               end do
           end if

           nam@$aName$ = rhs
           delete (rhs)
           delete (iComma)
       else
           nam@$aName$ = RHS    ; single string
       end if
       
       delete (iq)
       delete (CC)
       continue                 ; go to next iteration
   end if
                                ; MUST be integer or real
   delete (iq)                  ; iq referred to index of ' or "

   nDecimal = num(CC.eq.decimalc)      ; number of decimal pts

   if (nComma.gt.0) then
       iComma          = ind(CC.eq.commac)  ; inices of ,
       if (nDecimal.eq.0) then
           rhs         = new ( nComma+1, "integer")
           rhs(0)      = stringtointeger( chartostring( CC(iEq+1:iComma(0)-1) ))
           rhs(nComma) = stringtointeger( chartostring( CC(iComma(nComma-1)+1:nc-1) ))
           if (nComma.gt.1) then
               do i=0,nComma-2
                  rhs(i+1) = stringtointeger(chartostring(CC(iComma(i)+1:iComma(i+1)-1)))
               end do
           end if
       else
           rhs         = new ( nComma+1, "float")
           rhs(0)      = stringtofloat( chartostring( CC(iEq+1:iComma(0)-1) ))
           rhs(nComma) = stringtofloat( chartostring( CC(iComma(nComma-1)+1:nc-1) ))
           if (nComma.gt.1) then
               do i=0,nComma-2
                  rhs(i+1) = stringtofloat( chartostring(CC(iComma(i)+1:iComma(i+1)-1)))
               end do
           end if
       end if
       delete (rhs@_FillValue)
       delete (iComma)
   else
       if (nDecimal.eq.0) then
           rhs = stringtointeger(RHS)
       else
           rhs = stringtofloat(RHS)
       end if
   end if

   nam@$aName$ = rhs     ; associate values with variable
   delete (rhs)

   delete (CC)
  end do

  return(nam)
end





undef("keyword_values")
function keyword_values(pfName[1]:string, keynam[1]:string, return_type[1]:string)
;--------------------------------------------------
;  Read a 'simple'  ASCII (text): key=value(s)
;  Similar in structure to a Fortran NAMELIST
;--------------------------------------------------
; Nomenclature:
;    pfName - path for file: dir+filname
;    keynam - key name [no spaces; _ allowed)
;    delim  - delimiters/ seperators
;    return_type - "integer", "logical", "float", "string"
;
local types, tab, dq, delim, keyval, nline, KEY, iKEY, nKEY, nfld, chr, nchr, n, val, value
begin
;---Error test
   types = (/ "integer", "float", "double", "logical", "string" /)
   if (.not.(any(types.eq.return_type))) then 
       print("keyword_values:  return_type="+return_type+" not supported: ?correct spelling?")
       exit
   end if

;---Create delimiter string
   tab     = str_get_tab()        ; TAB
   dq      = str_get_dq()         ; "
   delim   = "=, "+dq+tab         ; add TAB and " 

;---read text (ascii) file and 'squeeze' leading spaces and TABs
   keyval = asciiread(pfName, -1, "string")        ; p(ath)/(f)ile Name
   keyval = str_squeeze(keyval)                    ; eliminate leading/trailing spaces and TABs
   nline  = dimsizes(keyval)                       ; # of lines (rows)

;---Extract field #1 (keynam)
   KEY    = str_get_field(keyval, 1, delim)        ; extract 'keynam' (leftmost) field 

;---Index of line containing 'keynam'
   iKEY   = ind(keynam.eq.KEY)                     ; index of line with 'keynam'
   if (ismissing(iKEY)) then
       print("keyword_values:  key name="+keynam+" not found: ?correct spelling?")
       return(default_fillvalue(return_type))
       exit
   end if

   nKEY   = dimsizes(iKEY)                          ; # keynam found             
   if (nKEY.gt.1) then
       print("keyword_values:  There are multiple entries named "+keynam)
       print("keyword_values:  The indices of the multiple entries are:")
       print("                 "+iKEY)
       return(default_fillvalue(return_type))
       exit
   end if

;---Number of fields for 'keyval(iKEY)' line.  This *includes* the keyword.
   nfld  = str_fields_count(keyval(iKEY), delim)   
   if (nfld.eq.1) then
       print("keyword_values:  key name="+keynam+" has no associated value(s)")
       return(default_fillvalue(return_type))
       exit
   end if

;---Create return field
   value = new(nfld-1, return_type)

   if (return_type.eq."string") then
       do n=1,nfld-1
          value(n-1) = str_get_field(keyval(iKEY), n+1, delim) 
       end do
   else
       if (return_type.eq."logical") then
           value = False                ; set default return value to False
           do n=1,nfld-1
              val = str_get_field(keyval(iKEY), n+1, delim)
              if (val.eq.".true." .or. val.eq.".TRUE." .or.  val.eq."TRUE" .or. \   
                  val.eq.".T." .or. val.eq.".t."       .or.  val.eq."T"    .or. \
                  val.eq."True".or. val.eq.".True.") then
                  value(n-1) = True     ; NCL logical
              end if   
           end do
        else
           do n=1,nfld-1
              value(n-1) = totype(str_get_field(keyval(iKEY), n+1, delim), return_type) 
           end do
        end if
   end if

   return(value) 
end
;----
