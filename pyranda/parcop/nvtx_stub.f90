! When nvtx is turned off, call these stubs instead 
module nvtx
   contains
   subroutine nvtxStartRange(tag, id)
      character(kind=c_char,len=*) :: tag
      integer, optional :: id
   end subroutine
   subroutine nvtxEndRange()
   end subroutine
end module nvtx
