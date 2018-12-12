pipeline {

    agent {
        docker {
            image 'infnpd/cmt-environment:latest-centos7'
            args '-u 0:0'
        }
    }

    stages {

        stage('Build') {

            steps {
                sh "yum -y localinstall http://artifacts.pd.infn.it/packages/MUOTOM/rpms/centos7/x86_64/base/cmt-ulib-0.2-1.el7.centos.x86_64.rpm"
                sh "yum -y localinstall http://artifacts.pd.infn.it/packages/MUOTOM/rpms/centos7/x86_64/base/cmt-ulib-devel-0.2-1.el7.centos.x86_64.rpm"
                sh "mkdir build && cd build && cmake -D ULIB_USE_QT5:BOOL=OFF .. && make"
            }

        }

    }

}
