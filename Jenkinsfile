pipeline {

    agent {
        docker {
            image 'infnpd/cmt-environment:latest-centos7'
            args '-u 0:0'
        }
    }

    parameters {
        string(name: 'repourl', defaultValue: 'http://artifacts.pd.infn.it/packages/MUOTOM/rpms/centos7/x86_64/base', description: 'Repository URL')
        string(name: 'ulibver', defaultValue: '0.2-1', description: 'uLib version')
    }

    stages {

        stage('Build') {

            steps {
                sh "yum -y localinstall ${params.repourl}/cmt-ulib-${params.ulibver}.el7.centos.x86_64.rpm"
                sh "yum -y localinstall ${params.repourl}/cmt-ulib-devel-${params.ulibver}.el7.centos.x86_64.rpm"
                sh "mkdir build && cd build && cmake -D ULIB_USE_QT5:BOOL=OFF .. && make"
            }

        }

    }

}
