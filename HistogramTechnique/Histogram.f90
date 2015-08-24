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
        double precision :: magnetizacao
        double precision :: magnetizacao2
        double precision :: magnetizacao4
        double precision :: energia
        double precision :: energia2
        double precision :: energia4
        double precision :: calorEspecifico
        double precision :: susceptibilidade
        double precision :: logMag2
        double precision :: cumulante
        double precision :: cumulanteE
    end type Medias
    !variaveis
    type(Medias) :: mediasTermodinamicas, somasTermodinamicas
    type(Configuracao) :: sistema

    double precision :: T0,T, J1, J2
    double precision :: delta !Delta Temperatura
    double precision :: eMax,eMin, argmax ! energias máximas e mínimas
    double precision :: energia,m   !energia e magnetizações
    double precision :: ex! = exp(n*(-delta*eSitio-argmax))

    double precision :: m2   = 0.0d0
    double precision :: m4   = 0.0d0
    double precision :: E    = 0.0d0
    double precision :: E2   = 0.0d0
    double precision :: E4   = 0.0d0

    double precision :: soma !soma dos ex
    double precision :: somaM  = 0.0d0
    double precision :: somaM2  = 0.0d0
    double precision :: somaM4  = 0.0d0

    double precision :: somaE    = 0.0d0
    double precision :: somaE2   = 0.0d0
    double precision :: somaE4   = 0.0d0

    double precision :: Mmed  =  0.0d0!somaMz/soma
    double precision :: M2med  =  0.0d0!somaMx2/soma
    double precision :: M4med  =  0.0d0!somaMx4/soma

    double precision :: Emed    =  0.0d0!somaE/soma
    double precision :: E2med   =  0.0d0!somaE2/soma
    double precision :: E4med   =  0.0d0!somaE2/soma

    integer :: contT, ContTF
    integer :: contAmostra
    integer :: n !numero de sitios
    integer :: numeroAmos
    integer :: i
    integer ::mag,  EJ1, EJ2

    open(1,file = 'dados.dat')
    open(3,file =  'magnetizacao.agr')
    open(4,file =  'magnetizacao2.agr')
    open(5,file =  'calorEspecifico.agr')
    open(6,file =  'susceptibilidade.agr')
    open(7,file =  'cumulante.agr')
    open(8,file =  'cumulanteE.agr')
    open(9,file =  'energia.agr')
    open(10,file =  'LogMag2.agr')

15  format (F18.7, 1x, F18.7, 1x, F18.7, 1x, F18.7)
111 format(F8.5,4X,E12.6)

    call lerDados(sistema)
    call gravaDados
    call maxmin
    J1=sistema%J1
    J2=sistema%J2
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
        somasTermodinamicas=medias(0.0d0, 0.0d00,  0.0d00, 0.0d00, 0.0d00, 0.0d00, 0.0d00, 0.0d00, 0.0d0, 0.0d00, 0.0d00)
        n=sistema%numeroSitios
        do contAmostra = 1, sistema%numeroAmostras
            soma    = 0.0d0
            somaM  = 0.0d0
            somaM2  = 0.0d0
            somaM4  = 0.0d0
            somaE    = 0.0d0
            somaE2   = 0.0d0
            somaE4   = 0.0d0
            do i=1, sistema%numeroPassosMC
                read(2,*)EJ1, EJ2, mag
                energia= J1*EJ1+ J2*EJ2
                ex      = exp(-delta*energia-argmax)
                m=mag/n
                m2     = m*m
                m4     = m2*m2
                E=energia
                E2=E*E
                E4=E2*E2
                soma   = soma   + ex
                somaM  = somaM + m*ex
                somaM2  = somaM2  + m2*ex
                somaM4  = somaM4  + m4*ex
                somaE    = somaE   + E*ex
                somaE2   = somaE2  + E2*ex
                somaE4   = somaE4  + E4*ex
            end do

            Mmed   = somaM/soma
            M2med  = somaM2/soma
            M4med  = somaM4/soma
            Emed    = somaE/soma
            E2med   = somaE2/soma
            E4med   = somaE4/soma

            somasTermodinamicas%magnetizacao     = somasTermodinamicas%magnetizacao   +  Mmed
            somasTermodinamicas%magnetizacao2   = somasTermodinamicas%magnetizacao2 +  M2med
            somasTermodinamicas%logMag2  = somasTermodinamicas%logMag2 +  log(M2med*n*n)
            somasTermodinamicas%energia           = somasTermodinamicas%energia         +  Emed
            somasTermodinamicas%calorEspecifico   = somasTermodinamicas%calorEspecifico +  (E2med - (Emed)**2)/(T*T)
            somasTermodinamicas%susceptibilidade = somasTermodinamicas%susceptibilidade + (M2med - (Mmed)**2)/T

            somasTermodinamicas%cumulante  = somasTermodinamicas%cumulante  +   (1-M4med/(3*M2med**2))
            somasTermodinamicas%cumulanteE  = somasTermodinamicas%cumulanteE  +   (1-E4med/(3*E2med**2))

        end do  !amostra
        close (2)
        numeroAmos=sistema%numeroAmostras
        mediasTermodinamicas%magnetizacao = somasTermodinamicas%magnetizacao/numeroAmos
        mediasTermodinamicas%magnetizacao2 = somasTermodinamicas%magnetizacao2/numeroAmos*n
        mediasTermodinamicas%logMag2 = somasTermodinamicas%logMag2/numeroAmos
        mediasTermodinamicas%energia = somasTermodinamicas%energia/numeroAmos
        mediasTermodinamicas%calorEspecifico = somasTermodinamicas%calorEspecifico/numeroAmos*n
        mediasTermodinamicas%susceptibilidade = somasTermodinamicas%susceptibilidade/numeroAmos*n
        mediasTermodinamicas%cumulante = somasTermodinamicas%cumulante/numeroAmos
        mediasTermodinamicas%cumulanteE = somasTermodinamicas%cumulanteE/numeroAmos

        write(*,*) T,mediasTermodinamicas%susceptibilidade
        write(3,111) T,mediasTermodinamicas%magnetizacao
        write(4,111) T,mediasTermodinamicas%magnetizacao2
        write(5,111) T,mediasTermodinamicas%calorEspecifico
        write(6,111) T,mediasTermodinamicas%susceptibilidade
        write(7,111) T,mediasTermodinamicas%cumulante
        write(8,111) T,mediasTermodinamicas%cumulanteE
        write(9,111) T,mediasTermodinamicas%energia
        write(10,111) T,mediasTermodinamicas%logMag2
        T = T + sistema%deltaT
    end do  !temperatura
!fim do programa

CONTAINS

    subroutine  maxmin
        !Variaveis Locais
        integer :: i_ , mcs_
        integer :: ej1_, ej2, m_
        double precision :: e_


15      format (F18.7, 1x, F18.7, 1x, F18.7, 1x, F18.7)
        open(2,file = 'hist.dat')
        read(2,*) ej1_, ej2, m_
        e_=J1*ej1+J2*ej2
        eMax = e_
        eMin = e_
        mcs_=sistema%numeroAmostras*sistema%numeroPassosMC
        do i_=1, mcs_-2
            read(2,15) ej1_, ej2, m_
            e_=J1*ej1+J2*ej2
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
        call get_command_argument(11, arg)
        read(arg,*) sistema_%tsimulacao
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
