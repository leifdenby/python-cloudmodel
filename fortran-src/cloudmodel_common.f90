module cloudmodel_common
   integer, parameter :: dp = selected_real_kind(12)
   integer :: nvars = 8

   ! Indexing into state-array
   integer :: i_r  = 1
   integer :: i_w  = 2
   integer :: i_T  = 3
   integer :: i_p  = 4
   integer :: i_qv = 5
   integer :: i_ql = 6
   integer :: i_qr = 7
   integer :: i_qi = 8

end module
