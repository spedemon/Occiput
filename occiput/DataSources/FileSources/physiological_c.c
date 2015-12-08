// occiput - 
// Stefano Pedemonte
// Harvard University, Boston. 
// July 2015


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <memory.h>

#define STATUS_SUCCESS              0
#define STATUS_IO_ERROR             1
#define STATUS_DECODE_ERROR         2 
#define STATUS_INITIALISATION_ERROR 3
#define STATUS_PARAMETER_ERROR      4
#define STATUS_UNHANDLED_ERROR      5 

#define BYTES_PER_PACKET            4
#define MAX_FILE_SIZE               10000000000   // 10 Gbytes
#define MAX_SAMPLES                 10000000      // 10 Msamples 


// Define the data structure 
#define PACKET_TYPE_DATA                0
#define PACKET_TYPE_CHANGE              1
#define PACKET_TYPE_TRIGGER_ECG         2
#define PACKET_TYPE_TRIGGER_BREATHING   4
#define PACKET_TYPE_TRIGGER_EXTERNAL    5
#define PACKET_TYPE_TRIGGER_PULSE_OX    6    
#define PACKET_TYPE_UNKNOWN             7

#define PACKET_SUBTYPE_UNKNOWN             100
#define PACKET_SUBTYPE_CHANGE_ECG_I        101
#define PACKET_SUBTYPE_CHANGE_ECG_aVF      102
#define PACKET_SUBTYPE_CHANGE_PULSE_OX     103
#define PACKET_SUBTYPE_CHANGE_BREATHING    104
#define PACKET_SUBTYPE_CHANGE_EXTERNAL     105
 
#define DATA_TYPE_UNKNOWN               9 
#define DATA_TYPE_ECG_I                 10
#define DATA_TYPE_ECG_aVF               11 
#define DATA_TYPE_PULSE_OX              12     
#define DATA_TYPE_BREATHING             13 
#define DATA_TYPE_EXTERNAL              14     

#define MASK_DATA                0b11111111 
#define MASK_TRIGGER_ECG         0b11111111 
#define MASK_TRIGGER_BREATHING   0b11111111 
#define MASK_TRIGGER_EXTERNAL    0b11111111
#define MASK_TRIGGER_PULSE_OX    0b11111111    
#define MASK_CHANGE              0b11111111 
#define MASK_CHANGE_ECG_I        0b11111111
#define MASK_CHANGE_ECG_aVF      0b11111111
#define MASK_CHANGE_PULSE_OX     0b11111111    
#define MASK_CHANGE_BREATHING    0b11111111 
#define MASK_CHANGE_EXTERNAL     0b11111111 

#define VALUE_DATA               0b00000000
#define VALUE_TRIGGER_ECG        0b01000000  
#define VALUE_TRIGGER_BREATHING  0b00010010  
#define VALUE_TRIGGER_EXTERNAL   0b00001000 
#define VALUE_TRIGGER_PULSE_OX   0b00100000     
#define VALUE_CHANGE             0b00000001
#define VALUE_CHANGE_ECG_I       0b00000001  
#define VALUE_CHANGE_ECG_aVF     0b00000010  
#define VALUE_CHANGE_PULSE_OX    0b00000101      
#define VALUE_CHANGE_BREATHING   0b00000110  
#define VALUE_CHANGE_EXTERNAL    0b00000111  
    
#define DATA_BITS                0x0000FFFF   //16 bits, right hand side 

    

/*
Return the status flags, so that they can be used by the (Python) wrapper to interpret 
the result of the function calls. 
*/
extern int status_success(int *status)
{
    *status = STATUS_SUCCESS; 
    return STATUS_SUCCESS; 
}
extern int status_io_error(int *status)
{
    *status = STATUS_IO_ERROR; 
    return STATUS_SUCCESS; 
}
extern int status_decode_error(int *status)
{
    *status = STATUS_DECODE_ERROR; 
    return STATUS_SUCCESS; 
}
extern int status_initialisation_error(int *status)
{
    *status = STATUS_INITIALISATION_ERROR; 
    return STATUS_SUCCESS; 
}
extern int status_parameter_error(int *status)
{
    *status = STATUS_PARAMETER_ERROR; 
    return STATUS_SUCCESS; 
}
extern int status_unhandled_error(int *status)
{
    *status = STATUS_UNHANDLED_ERROR; 
    return STATUS_SUCCESS; 
}


/*
Return the integer parameter. This function is meant to test the (Python) wrapper of the library. 
*/
extern int echo(int *in, int *out)
{
    *out = *in; 
    return STATUS_SUCCESS; 
}




/* 
Print n as a binary number. Utility function for debugging. 
*/
void printbits(int n) 
{
    unsigned int i, step;

    if (0 == n)  /* For simplicity's sake, I treat 0 as a special case*/
    {
        printf("0000");
        return;
    }

    i = 1<<(sizeof(n) * 8 - 1);

    step = -1; /* Only print the relevant digits */
    step >>= 4; /* In groups of 4 */
    while (step >= n) 
    {
        i >>= 4;
        step >>= 4;
    }

    /* At this point, i is the smallest power of two larger or equal to n */
    while (i > 0) 
    {
        if (n & i)
            printf("1");
        else
            printf("0");
        i >>= 1;
    }
    printf("\n");
}


typedef struct Packet {
   char buffer[4];  
   int type; 
   int subtype; 
   int data; 
} Packet; 


int clear_packet(Packet *packet)
{
   packet->type            = PACKET_TYPE_UNKNOWN; 
   packet->subtype         = PACKET_SUBTYPE_UNKNOWN; 
   packet->data            = 0;
   return STATUS_SUCCESS;
}



int decode_packet(Packet* packet)
{
    // Uncomment the following lines to print packet information for debugging. 
    //int *int_p; 
    //int_p = (int*) packet->buffer; 
    //printbits(*int_p); 

    //###### Detect the type of packet ######
    //detect bits configuration: check if the 'binary and' of PACKET and MASK equals VALUE: 
    if ( !((packet->buffer[3] & MASK_DATA) - VALUE_DATA) ) {
        packet->type = PACKET_TYPE_DATA;  
        packet->subtype = PACKET_SUBTYPE_UNKNOWN;
    }
    else if ( !((packet->buffer[3] & MASK_TRIGGER_ECG) - VALUE_TRIGGER_ECG) ) {
        packet->type = PACKET_TYPE_TRIGGER_ECG;  
        packet->subtype = PACKET_SUBTYPE_UNKNOWN;
    }
    else if ( !((packet->buffer[3] & MASK_TRIGGER_BREATHING) - VALUE_TRIGGER_BREATHING) ) {
        packet->type = PACKET_TYPE_TRIGGER_BREATHING;  
        packet->subtype = PACKET_SUBTYPE_UNKNOWN;
    }
    else if ( !((packet->buffer[3] & MASK_TRIGGER_EXTERNAL) - VALUE_TRIGGER_EXTERNAL) ) {
        packet->type = PACKET_TYPE_TRIGGER_EXTERNAL; 
        packet->subtype = PACKET_SUBTYPE_UNKNOWN;
    }
    else if ( !((packet->buffer[3] & MASK_TRIGGER_PULSE_OX) - VALUE_TRIGGER_PULSE_OX) ) {
        packet->type = PACKET_TYPE_TRIGGER_PULSE_OX; 
        packet->subtype = PACKET_SUBTYPE_UNKNOWN;
    }
    else if ( !((packet->buffer[3] & MASK_CHANGE) - VALUE_CHANGE) ) {
        packet->type = PACKET_TYPE_CHANGE; 
        // decode subtype
        if ( !((packet->buffer[2] & MASK_CHANGE_ECG_I) - VALUE_CHANGE_ECG_I) ) {
            packet->subtype = PACKET_SUBTYPE_CHANGE_ECG_I;
        }
        else if ( !((packet->buffer[2] & MASK_CHANGE_ECG_aVF) - VALUE_CHANGE_ECG_aVF) ) {
            packet->subtype = PACKET_SUBTYPE_CHANGE_ECG_aVF;
        }
        else if ( !((packet->buffer[2] & MASK_CHANGE_PULSE_OX) - VALUE_CHANGE_PULSE_OX) ) {
            packet->subtype = PACKET_SUBTYPE_CHANGE_PULSE_OX;
        }
        else if ( !((packet->buffer[2] & MASK_CHANGE_BREATHING) - VALUE_CHANGE_BREATHING) ) {
            packet->subtype = PACKET_SUBTYPE_CHANGE_BREATHING; 
        }
        else if ( !((packet->buffer[2] & MASK_CHANGE_EXTERNAL) - VALUE_CHANGE_EXTERNAL) ) {
            packet->subtype = PACKET_SUBTYPE_CHANGE_EXTERNAL;
        }
    }
    else {
        packet->type = PACKET_TYPE_UNKNOWN; 
        packet->subtype = PACKET_SUBTYPE_UNKNOWN;
    }

    //###### Read the content of the packet ######
    //## Only data packets have data (is it true? FIXME) ##
    //if (packet->type == PACKET_TYPE_DATA) {
    packet->data = ( *(int*)(packet->buffer) & DATA_BITS );  
    //    }

   return STATUS_SUCCESS; 
}


int read_packets(FILE *stream,char *buffer, int64_t number_of_packets)
{
    if( fread(buffer, BYTES_PER_PACKET, number_of_packets, stream ) != number_of_packets)   
        return STATUS_IO_ERROR;  
    return STATUS_SUCCESS; 
}


int read_packet(FILE *stream, Packet* packet) 
{
    if( fread(packet->buffer,4,1,stream) != 1) 
        return STATUS_IO_ERROR; 
    return STATUS_SUCCESS; 
}


int read_and_decode_packet(FILE *stream, Packet* packet)
{
    if( read_packet(stream, packet) != STATUS_SUCCESS) 
        return STATUS_IO_ERROR; 
    return decode_packet(packet);
}


unsigned short G_samples_ecg_I[MAX_SAMPLES]; 
unsigned short G_samples_ecg_aVF[MAX_SAMPLES]; 
unsigned short G_samples_pulse_ox[MAX_SAMPLES]; 
unsigned short G_samples_breathing[MAX_SAMPLES]; 
unsigned short G_samples_external[MAX_SAMPLES]; 


unsigned short G_triggers_ecg_I[MAX_SAMPLES]; 
unsigned short G_triggers_ecg_aVF[MAX_SAMPLES];
unsigned short G_triggers_pulse_ox[MAX_SAMPLES];
unsigned short G_triggers_breathing[MAX_SAMPLES]; 
unsigned short G_triggers_external[MAX_SAMPLES]; 

//unsigned short G_unknown_packets[MAX_SAMPLES]; 

unsigned int G_n_samples_ecg_I     = 0;
unsigned int G_n_samples_ecg_aVF   = 0;
unsigned int G_n_samples_pulse_ox  = 0;
unsigned int G_n_samples_breathing = 0;
unsigned int G_n_samples_external  = 0;
//unsigned int G_n_unknown           = 0;



void clear_global_data()
{
    G_n_samples_ecg_I        = 0;
    G_n_samples_ecg_aVF      = 0;
    G_n_samples_pulse_ox     = 0;
    G_n_samples_breathing    = 0;
    G_n_samples_external     = 0; 
//    G_n_unknown              = 0;
    memset((void *) G_samples_ecg_I, 0, sizeof(unsigned short)*MAX_SAMPLES );
    memset((void *) G_samples_ecg_aVF, 0, sizeof(unsigned short)*MAX_SAMPLES );
    memset((void *) G_samples_pulse_ox, 0, sizeof(unsigned short)*MAX_SAMPLES );
    memset((void *) G_samples_breathing, 0, sizeof(unsigned short)*MAX_SAMPLES );
    memset((void *) G_samples_external, 0, sizeof(unsigned short)*MAX_SAMPLES );
    
    memset((void *) G_triggers_ecg_I, 0, sizeof(unsigned short)*MAX_SAMPLES );
    memset((void *) G_triggers_ecg_aVF, 0, sizeof(unsigned short)*MAX_SAMPLES );
    memset((void *) G_triggers_pulse_ox, 0, sizeof(unsigned short)*MAX_SAMPLES );
    memset((void *) G_triggers_breathing, 0, sizeof(unsigned short)*MAX_SAMPLES );
    memset((void *) G_triggers_external, 0, sizeof(unsigned short)*MAX_SAMPLES );    
//    memset((void *) G_unknown_packets, 0, sizeof(unsigned short)*MAX_SAMPLES );  
}




// Functions exported to (Python) wrapper start here. //

extern int load_data(char *filename)
{
    clear_global_data(); 
    
    unsigned int status = STATUS_SUCCESS; 
    unsigned int i,j; 
    unsigned int done = 0; 
    unsigned int buffer_size_batch = 1; //change this to load batches of packets for increased speed 
    unsigned int n_packets_batch   = 1; //change this to load batches of packets for increased speed
    char batch_buffer[buffer_size_batch*BYTES_PER_PACKET];
    Packet packet; 
    
    unsigned int current_data_type = DATA_TYPE_UNKNOWN;  //this changes when type-change packets are received 

    FILE *fid; 
    fid=fopen(filename, "rb");
    if (fid == NULL) {
        fprintf(stderr,"Failed to open binary physio file. \n");
        status = STATUS_IO_ERROR; 
        return status; 
    }
    for (j=0; j<MAX_FILE_SIZE; j++) {
        status = read_packets(fid, batch_buffer, n_packets_batch); 
        if (status!=STATUS_SUCCESS) {
            fprintf(stderr,"Problem reading the physio binary data (batch index: %d  batch size: %d).\n",j,n_packets_batch); 
            //return status; 
            status = STATUS_SUCCESS; 
            done = 1; 
            break; 
        }
        for (i = 0; i < n_packets_batch; i++) { 
            clear_packet(&packet);
            memcpy ( packet.buffer, batch_buffer+i*BYTES_PER_PACKET, BYTES_PER_PACKET );
            status = decode_packet(&packet);
            if (status!=STATUS_SUCCESS) {
                fprintf(stderr,"Problem decoding the listmode data.\n"); 
                fclose(fid); 
                return status; 
            }    
            if (packet.type == PACKET_TYPE_DATA) {
                if (current_data_type==DATA_TYPE_ECG_I) {
                    G_samples_ecg_I[G_n_samples_ecg_I]=packet.data; 
                    G_n_samples_ecg_I+=1; 
                }
                else if (current_data_type==DATA_TYPE_ECG_aVF) {
                    G_samples_ecg_aVF[G_n_samples_ecg_aVF]=packet.data; 
                    G_n_samples_ecg_aVF+=1; 
                }
                else if (current_data_type==DATA_TYPE_PULSE_OX) {
                    G_samples_pulse_ox[G_n_samples_pulse_ox]=packet.data; 
                    G_n_samples_pulse_ox+=1; 
                }
                else if (current_data_type==DATA_TYPE_BREATHING) {
                    G_samples_breathing[G_n_samples_breathing]=packet.data;
                    G_n_samples_breathing+=1;
                }
                else if (current_data_type==DATA_TYPE_EXTERNAL) {
                    G_samples_external[G_n_samples_external]=packet.data;
                    G_n_samples_external+=1;
                }
                else {
                    //G_n_samples_unknown+=1; 
                }
            }
            else if (packet.type == PACKET_TYPE_CHANGE) { 
//                fprintf(stderr,"type change ecg\n");
                if (packet.subtype == PACKET_SUBTYPE_CHANGE_ECG_I)
                    current_data_type = DATA_TYPE_ECG_I; 
                else if (packet.subtype == PACKET_SUBTYPE_CHANGE_ECG_aVF)
                    current_data_type = DATA_TYPE_ECG_aVF; 
                else if (packet.subtype == PACKET_SUBTYPE_CHANGE_PULSE_OX)
                    current_data_type = DATA_TYPE_PULSE_OX;                 
                else if (packet.subtype == PACKET_SUBTYPE_CHANGE_BREATHING)
                    current_data_type = DATA_TYPE_BREATHING; 
                else if (packet.subtype == PACKET_SUBTYPE_CHANGE_EXTERNAL)
                    current_data_type = DATA_TYPE_EXTERNAL; 
                else {
                }
            }
            else if (packet.type == PACKET_TYPE_TRIGGER_ECG) { 
                //fprintf(stderr,"type trigger ecg\n");
                if (current_data_type == DATA_TYPE_ECG_I) { 
                    G_triggers_ecg_I[G_n_samples_ecg_I]=packet.data;
                    G_samples_ecg_I[G_n_samples_ecg_I]=packet.data;
                    G_n_samples_ecg_I+=1;
                }
                else if (current_data_type == DATA_TYPE_ECG_aVF) { 
                    G_triggers_ecg_aVF[G_n_samples_ecg_aVF]=packet.data;
                    G_samples_ecg_aVF[G_n_samples_ecg_aVF]=packet.data;
                    G_n_samples_ecg_aVF+=1;
                }                
            }
            else if (packet.type == PACKET_TYPE_TRIGGER_BREATHING) { 
                //fprintf(stderr,"type trigger breath\n");
                G_triggers_breathing[G_n_samples_breathing]=packet.data;
                G_samples_breathing[G_n_samples_breathing]=packet.data;
                G_n_samples_breathing+=1;
            }
            else if (packet.type == PACKET_TYPE_TRIGGER_EXTERNAL) { 
                //fprintf(stderr,"type trigger ext\n");
                G_triggers_external[G_n_samples_external]=packet.data;
                G_samples_external[G_n_samples_external]=packet.data;
                G_n_samples_external+=1; 
            }
            else if (packet.type == PACKET_TYPE_TRIGGER_PULSE_OX) { 
                //fprintf(stderr,"type trigger ext\n");
                G_triggers_pulse_ox[G_n_samples_pulse_ox]=packet.data;
                G_samples_pulse_ox[G_n_samples_pulse_ox]=packet.data;
                G_n_samples_pulse_ox+=1;
            }
            else {
                //G_unknown_packets[]=packet.data;
            }
          
        }
    }
    return status; 
}



                
extern int get_info(unsigned int *n_samples_breathing, unsigned int *n_samples_ecg_I, unsigned int *n_samples_ecg_aVF, unsigned int *n_samples_pulse_ox, unsigned int *n_samples_external)
{
        *n_samples_breathing = G_n_samples_breathing;
        *n_samples_ecg_I     = G_n_samples_ecg_I;
        *n_samples_ecg_aVF   = G_n_samples_ecg_aVF;
        *n_samples_pulse_ox  = G_n_samples_pulse_ox;
        *n_samples_external  = G_n_samples_external;
        return 0; 
}


extern int get_data(
    unsigned int *n_samples_breathing, unsigned int *n_samples_ecg_I, unsigned int *n_samples_ecg_aVF, 
    unsigned int *n_samples_pulse_ox, unsigned int *n_samples_external, 
    unsigned short *samples_breathing, unsigned short *samples_ecg_I, unsigned short *samples_ecg_aVF, 
    unsigned short *samples_pulse_ox, unsigned short *samples_external, 
    unsigned short *triggers_breathing, unsigned short *triggers_ecg_I, unsigned short *triggers_ecg_aVF, 
    unsigned short *triggers_pulse_ox, unsigned short *triggers_external)
{
    unsigned int N_breathing = *n_samples_breathing; 
    unsigned int N_ecg_I     = *n_samples_ecg_I; 
    unsigned int N_ecg_aVF   = *n_samples_ecg_aVF; 
    unsigned int N_pulse_ox  = *n_samples_pulse_ox; 
    unsigned int N_external  = *n_samples_external; 
    
    memcpy ( (void*) samples_breathing,  (void*) G_samples_breathing,    sizeof(unsigned short)*N_breathing );
    memcpy ( (void*) samples_ecg_I,      (void*) G_samples_ecg_I,        sizeof(unsigned short)*N_ecg_I );
    memcpy ( (void*) samples_ecg_aVF,    (void*) G_samples_ecg_aVF,      sizeof(unsigned short)*N_ecg_aVF );
    memcpy ( (void*) samples_pulse_ox,   (void*) G_samples_pulse_ox,     sizeof(unsigned short)*N_pulse_ox );
    memcpy ( (void*) samples_external,   (void*) G_samples_external,     sizeof(unsigned short)*N_external );
    memcpy ( (void*) triggers_breathing, (void*) G_triggers_breathing,   sizeof(unsigned short)*N_breathing );
    memcpy ( (void*) triggers_ecg_I,     (void*) G_triggers_ecg_I,       sizeof(unsigned short)*N_ecg_I );
    memcpy ( (void*) triggers_ecg_aVF,   (void*) G_triggers_ecg_aVF,     sizeof(unsigned short)*N_ecg_aVF );
    memcpy ( (void*) triggers_pulse_ox,  (void*) G_triggers_pulse_ox,    sizeof(unsigned short)*N_pulse_ox );
    memcpy ( (void*) triggers_external,  (void*) G_triggers_external,    sizeof(unsigned short)*N_external );    
    return 0; 
}






