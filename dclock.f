      function dclock
      implicit none
      real*8::dclock
      integer::clk,clk_rate,clk_max

      call system_clock(clk,clk_rate,clk_max)
      dclock = real(clk,8)/clk_rate
      end function dclock
