program Histogram
    implicit none

    type :: Configuracao
        integer :: l
        integer :: numeroCoordenacao
        integer :: numeroSitios
        integer :: numeroTermalizacao
        integer :: numeroPassosMC
        integer :: numeroAmostras
        integer :: numeroIntervalosTeperatura
        double precision :: tInicial
        double precision :: tFinal
        double precision :: deltaT
        double precision :: concentracao
        double precision :: tSimulacao
        double precision :: J1
        double precision :: J2
    end type Configuracao

    type :: Medias
        double precision :: magnetizacaoZ
        double precision :: magnetizacaoX
        double precision :: magnetizacaoX2
        double precision :: magnetizacaoXY
        double precision :: magnetizacaoXY2
        double precision :: magnetizacaoXY4
        double precision :: energia
        double precision :: energia2
        double precision :: energia4
        double precision :: helicidade1
        double precision :: helicidade2
        double precision :: helicidade3
        double precision :: calorEspecifico
        double precision :: susceptibilidadeX
        double precision :: susceptibilidadeXY
        double precision :: moduloHelicidade
        double precision :: logMagXY2
        double precision :: cumulantex
        double precision :: cumulantexy
        double precision :: cumulanteE
    end type Medias
    !variaveis
    type(Medias) :: mediasTermodinamicas, somasTermodinamicas
    type(Configuracao) :: sistema

    double precision :: T0,T
    double precision :: delta !Delta Temperatura
    double precision :: eMax,eMin, argmax ! energias máximas e mínimas
    double precision :: eSitio,mx,my,mz   !energia e magnetizações
    double precision :: ex! = exp(n*(-delta*eSitio-argmax))

    double precision :: mx2!     = mx*mx
    double precision :: mx4!     = mx2*mx2
    double precision :: mxy2!    = mx2+my*my
    double precision :: mxy4!    = mxy2*mxy2
    double precision :: mxy!     = sqrt(mxy2)
    double precision :: E    = 0.0d0
    double precision :: E2   = 0.0d0
    double precision :: E4   = 0.0d0

    double precision :: soma !soma dos ex
    double precision :: somaMz  = 0.0d0
    double precision :: somaMx   = 0.0d0
    double precision :: somaMxy  = 0.0d0
    double precision :: somaMx2  = 0.0d0
    double precision :: somaMxy2 = 0.0d0
    double precision :: somaMx4  = 0.0d0
    double precision :: somaMxy4 = 0.0d0
    double precision :: somaE    = 0.0d0
    double precision :: somaE2   = 0.0d0
    double precision :: somaE4   = 0.0d0

    double precision :: Mzmed  =  0.0d0!somaMz/soma
    double precision :: Mxmed   =  0.0d0!somaMx/soma
    double precision :: Mxymed  = 0.0d0! somaMxy/soma
    double precision :: Mx2med  =  0.0d0!somaMx2/soma
    double precision :: Mxy2med =  0.0d0!somaMxy2/soma
    double precision :: Mx4med  =  0.0d0!somaMx4/soma
    double precision :: Mxy4med =  0.0d0!somaMxy4/soma
    double precision :: Emed    =  0.0d0!somaE/soma
    double precision :: E2med   =  0.0d0!somaE2/soma
    double precision :: E4med   =  0.0d0!somaE2/soma

    integer :: contT, ContTF
    integer :: contAmostra
    integer :: n !numero de sitios
    integer :: numeroAmos
    integer :: i


    open(1,file = 'dados.dat')

    open(3,file = 'magnetizacaoZ.agr')
    open(4,file = 'magnetizacaoXY.agr')
    open(5,file = 'magnetizacaoXY2.agr')

    open(7,file = 'calorEspecifico.agr')
    open(8,file = 'susceptibilidadeX.agr')
    open(9,file = 'susceptibilidadeXY.agr')
    open(10,file = 'cumulantex.agr')
    open(11,file = 'cumulantexy.agr')
    open(12,file = 'cumulanteE.agr')
    open(13,file = 'energia.agr')
    open(14,file = 'magnetizacaox.agr')
    open(14,file = 'LogMagXY2.agr')

15  format (F18.7, 1x, F18.7, 1x, F18.7, 1x, F18.7)
111 format(F8.5,4X,E12.6)
    call lerDados(sistema)
    call gravaDados
    call maxmin
    T  = sistema%tInicial
    T0 = sistema%tsimulacao
    n  = sistema%numeroSitios
    ContTF = sistema%numeroIntervalosTeperatura
    write(*,*) "ok"
    do contT=1,ContTF  !repetição na temperatura
        delta = 1.0/T -1.0/T0

        if (delta > 0) then
            argmax = -delta*(eMin)
        else
            argmax = -delta*(eMax)
        end if
        open(2,file = 'hist.dat')
        somasTermodinamicas=medias(0.0d0, 0.0d00,  0.0d00, 0.0d00, 0.0d00, 0.0d00, 0.0d00, 0.0d00,0.0d0, 0.0d00, &
            &0.0d00, 0.0d00, 0.0d00, 0.0d00,0.0d0,  0.0d00, 0.0d00, 0.0d00, 0.0d00,0.0d00)
        n=sistema%numeroSitios
        do contAmostra = 1, sistema%numeroAmostras
            soma    = 0.0d0
            somaMz  = 0.0d0
            somaMx   = 0.0d0
            somaMxy  = 0.0d0
            somaMx2  = 0.0d0
            somaMxy2 = 0.0d0
            somaMx4  = 0.0d0
            somaMxy4 = 0.0d0
            somaE    = 0.0d0
            somaE2   = 0.0d0
            somaE4   = 0.0d0
            do i=1, sistema%numeroPassosMC
                read(2,15) m, EJ1, EJ2
                energia=j1*EJ1+ J2*EJ2
                m/n
                ex      = exp((-delta*eSitio-argmax))
                mx2     = mx**2
                mx4     = mx2**2
                mxy2    = mx2+my**2
                mxy4    = mxy2**2
                mxy     = sqrt(mxy2)
                E=eSitio*n
                E2=E**2
                E4=E2**2
                soma   = soma   + ex
                somaMz  = somaMz + mz*ex
                somaMx   = somaMx + mx*ex
                somaMxy  = somaMxy + mxy*ex
                somaMx2  = somaMx2  + mx2*ex
                somaMxy2 = somaMxy2 + mxy2*ex
                somaMx4  = somaMx4  + mx4*ex
                somaMxy4 = somaMxy4 + mxy4*ex
                somaE    = somaE   + E*ex
                somaE2   = somaE2  + E2*ex
                somaE2   = somaE4  + E4*ex
            end do

            Mzmed  = somaMz/soma
            Mxmed   = somaMx/soma
            Mxymed  = somaMxy/soma
            Mx2med  = somaMx2/soma
            Mxy2med = somaMxy2/soma
            Mx4med  = somaMx4/soma
            Mxy4med = somaMxy4/soma
            Emed    = somaE/soma
            E2med   = somaE2/soma
            E4med   = somaE4/soma

            somasTermodinamicas%magnetizacaoZ     = somasTermodinamicas%magnetizacaoZ   +  MzMed
            somasTermodinamicas%magnetizacaoXY    = somasTermodinamicas%magnetizacaoXY  +  Mxymed
            somasTermodinamicas%magnetizacaoX    = somasTermodinamicas%magnetizacaoX  +  Mxmed

            somasTermodinamicas%magnetizacaoXY2   = somasTermodinamicas%magnetizacaoXY2 +  Mxy2med
            somasTermodinamicas%logMagXY2  = somasTermodinamicas%logMagXY2 +  log(Mxy2med*n*n)
            somasTermodinamicas%energia           = somasTermodinamicas%energia         +  Emed
            somasTermodinamicas%calorEspecifico   = somasTermodinamicas%calorEspecifico +  (E2med - (Emed)**2)/(T*T)
            somasTermodinamicas%susceptibilidadeX = somasTermodinamicas%susceptibilidadeX + (Mx2med - (Mxmed)**2)/T
            somasTermodinamicas%susceptibilidadeXY = somasTermodinamicas%susceptibilidadeXY + (Mxy2med - (Mxymed)**2)/T
            somasTermodinamicas%cumulantex  = somasTermodinamicas%cumulanteX  +   (1-Mx4med/(3*Mx2med**2))
            somasTermodinamicas%cumulantexy = somasTermodinamicas%cumulantexy +   (1-Mxy4med/(3*Mxy2med**2))
            somasTermodinamicas%cumulanteE  = somasTermodinamicas%cumulanteE  +   (1-E4med/(3*E2med**2))

        end do  !amostra
        close (2)
        numeroAmos=sistema%numeroAmostras
        mediasTermodinamicas%magnetizacaoZ = somasTermodinamicas%magnetizacaoZ/numeroAmos
        mediasTermodinamicas%magnetizacaox = somasTermodinamicas%magnetizacaox/numeroAmos
        mediasTermodinamicas%magnetizacaoXY = somasTermodinamicas%magnetizacaoXY/numeroAmos
        mediasTermodinamicas%magnetizacaoXY2 = somasTermodinamicas%magnetizacaoXY2/numeroAmos*n
        mediasTermodinamicas%logMagXY2 = somasTermodinamicas%logMagXY2/numeroAmos
        mediasTermodinamicas%energia = somasTermodinamicas%energia/numeroAmos
        mediasTermodinamicas%calorEspecifico = somasTermodinamicas%calorEspecifico/numeroAmos*n
        mediasTermodinamicas%susceptibilidadeX = somasTermodinamicas%susceptibilidadeX/numeroAmos*n
        mediasTermodinamicas%susceptibilidadeXY = somasTermodinamicas%susceptibilidadeXY/numeroAmos*n
        mediasTermodinamicas%cumulantex = somasTermodinamicas%cumulanteX/numeroAmos
        mediasTermodinamicas%cumulantexy = somasTermodinamicas%cumulantexy/numeroAmos
        mediasTermodinamicas%cumulanteE = somasTermodinamicas%cumulanteE/numeroAmos

        write(*,*) T,mediasTermodinamicas%susceptibilidadeXY
        write(3,111) T,mediasTermodinamicas%magnetizacaoZ
        write(4,111) T,mediasTermodinamicas%magnetizacaoXY
        write(5,111) T,mediasTermodinamicas%magnetizacaoXY2
        write(7,111) T,mediasTermodinamicas%calorEspecifico
        write(8,111) T,mediasTermodinamicas%susceptibilidadeX
        write(9,111)  T,mediasTermodinamicas%susceptibilidadeXY
        write(10,111) T,mediasTermodinamicas%cumulantex
        write(11,111) T,mediasTermodinamicas%cumulantexy
        write(12,111) T,mediasTermodinamicas%cumulanteE
        write(13,111) T,mediasTermodinamicas%energia
        write(14,111) T,mediasTermodinamicas%magnetizacaoX
        write(15,111) T,mediasTermodinamicas%logMagXY2

        T = T + sistema%deltaT
    end do  !temperatura
!fim do programa

CONTAINS
    subroutine  maxmin
        !Variaveis Locais
        integer :: i_ , mcs_
        double precision :: e_, mx_, my_,mz_

15      format (F18.7, 1x, F18.7, 1x, F18.7, 1x, F18.7)
        open(2,file = 'hist.dat')
        read(2,15) e_,mx_,my_,mz_
        eMax = e_
        eMin = e_
        mcs_=sistema%numeroAmostras*sistema%numeroPassosMC
        do i_=1, mcs_-2
            read(2,15) e_,mx_,my_,mz_
            if (e_ > eMax) eMax = e_
            if (e_ < eMin) eMin = e_
        end do
        close(2)
        return
    end subroutine
    !-----------------------------------------------------------------------------

    subroutine gravaDados
        !Variaveis Locais
        integer :: i_

        do i_ = 3, 12
            write(i_,*) "# L : ", sistema%l
            write(i_,*) "# numero de coordenação : ", sistema%numeroCoordenacao
            write(i_,*) "# Numero de passos para termalização : ", sistema%numeroTermalizacao
            write(i_,*) "# Número de passos MC : ", sistema%numeroPassosMC
            write(i_,*) "# Concentração : ", sistema%concentracao
            write(i_,*) "# Temperatura inicial : ", sistema%tInicial
            write(i_,*) "# Temperatura final : ", sistema%tFinal
            write(i_,*) "# delta T : ", sistema%deltaT
        end do

    end subroutine gravaDados
    !-----------------------------------------------------------------------------
    subroutine leargumentos(sistema_)
        type(configuracao), intent(out) :: sistema_
        character(len=32) :: arg
        call get_command_argument(1, arg)
        read (arg,*) sistema_%l
        call get_command_argument(2, arg)
        read(arg,*) sistema_%concentracao
        call get_command_argument(3, arg)
        read(arg,*) sistema_%tInicial
        call get_command_argument(4, arg)
        read(arg,*) sistema_%tFinal
        call get_command_argument(5, arg)
        read(arg,*) sistema_%deltaT
        call get_command_argument(6, arg)
        read(arg,*) sistema_%J1
        call get_command_argument(7, arg)
        read(arg,*) sistema_%J2
        call get_command_argument(8, arg)
        read(arg,*) sistema_%numeroTermalizacao
        call get_command_argument(9, arg)
        read(arg,*) sistema_%numeroPassosMC
        call get_command_argument(10, arg)
        read(arg,*) sistema_%numeroAmostras
        sistema_%numeroSitios = 2*sistema_%L*sistema_%L*sistema_%L
        sistema_%numeroIntervalosTeperatura = int((sistema_%tFinal -sistema_%tInicial)/sistema_%deltaT) +1
    end subroutine leargumentos
      !-----------------------------------------------------------------------------
    subroutine lerDados(sistema_)
        !Variaveis mudas
        type(configuracao), intent(out) :: sistema_

        !Variaveis Locais
        character*20 :: lendo

        read(1,*)  lendo
        read(1,*)  sistema_%l
        write(*,*) lendo, sistema_%l

        read(1,*)  lendo
        read(1,*)  sistema_%numeroCoordenacao
        write(*,*) lendo, sistema_%numeroCoordenacao

        read(1,*)  lendo
        read(1,*)  sistema_%numeroTermalizacao
        write(*,*) lendo, sistema_%numeroTermalizacao

        read(1,*)  lendo
        read(1,*)  sistema_%numeroPassosMC
        write(*,*) lendo, sistema_%numeroPassosMC

        read(1,*)  lendo
        read(1,*)  sistema_%concentracao
        write(*,*) lendo, sistema_%concentracao

        read(1,*)  lendo
        read(1,*)  sistema_%tInicial
        write(*,*) lendo, sistema_%tInicial

        read(1,*)  lendo
        read(1,*)  sistema_%tFinal
        write(*,*) lendo, sistema_%tFinal

        read(1,*)  lendo
        read(1,*)  sistema_%deltaT
        write(*,*) lendo, sistema_%deltaT

        read(1,*)  lendo
        read(1,*)  sistema_%numeroAmostras
        write(*,*) lendo, sistema_%numeroAmostras


        read(1,*)  lendo
        read(1,*)  sistema_%tsimulacao
        write(*,*) lendo, sistema_%tsimulacao


        sistema_%numeroSitios =2*sistema_%L*sistema_%L*sistema_%L
        sistema_%numeroIntervalosTeperatura = int((sistema_%tFinal -sistema_%tInicial)/sistema_%deltaT) +1

    end subroutine LerDados

end program Histogram
