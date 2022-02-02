#include "icgemio.hpp"
#include <cstdio>

int main(int argc, char *argv[]) {
    if (argc!=2 && argc!=3) {
        fprintf(stderr, "Usage: %s <icgem-file> [MAX DEGREE (optional)]\n", argv[0]);
        return 1;
    }

    // an icgem (gfc) instance
    Icgem gfc(argv[1]);
    
    // parse the header ...
    if (gfc.parse_header()) {
        fprintf(stderr, "ERROR! Failed to parse icgem header!\n");
        return 1;
    }

    // print model info
    #ifdef DEBUG 
    else {
        gfc.print_details();
    }
    #endif

    int degree = gfc.degree();
    if (argc==3) degree = std::atoi(argv[2]);

    // allocate memory to store harmonic coefficients
    HarmonicCoeffs hc(degree);

    // parse data
    if (gfc.parse_data(degree, degree, &hc)) {
        fprintf(stderr, "ERROR! Failed to parsed harmonic coefficients\n");
        return 1;
    }

    char yn;
    // print coefficients
    #ifdef DEBUG
    printf("Want to see the coefficients [y/n]?\n");
    std::scanf("%c", &yn);
    if (yn=='y' || yn=='Y') hc.print();
    #endif

    bool ask_again = true;
    while (ask_again) {
        printf("Enter degree order to print Cnm and Snm values; if degree and/or order is negative program exits:");
        int dg, rd;
        std::scanf("%d %d", &dg, &rd);
        if (dg<0 || rd<0) break;
        if (dg>degree || rd >dg) {
            fprintf(stderr, "Erronuous degree and/or order value!");
        } else {
            printf("C(%d,%d)=%20.15e S(%d,%d)=%20.15e\n", dg,rd,hc.C(dg,rd),dg,rd,hc.S(dg,rd));
        }
    }

    return 0;
}