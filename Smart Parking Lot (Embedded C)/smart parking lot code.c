#include<reg51.h>
#include<string.h>
//password
sbit pw1=0x80;
sbit pw2=0x81;
sbit pw3=0x82;
sbit pw4=0x83;
//LCD
sbit rs = P3^4;  
sbit rw = P3^5;  
sbit en = P3^6;  
sfr data_lcd=0x90;
//motor for door
sfr doors=0xA0;
unsigned int slots=5;
unsigned char msg[]="PARKING FULL";
unsigned char str1[] = "Number of slots";
unsigned char str2[] = "available: ";
unsigned char pwmsg[] = "Enter Password";
unsigned char wrong[] = "Wrong Password";
void lcd_init();
void cmd_sr(unsigned char);
void data_sr(unsigned char);
void delay(){
	unsigned int i;
	for(i=0;i<10000;i++);
}

void cmd_sr(unsigned char x){
	data_lcd=x;
	rw=0;
	rs=0;
	en=1;
	delay();
	en=0
}

void data_sr(unsigned char x){
	data_lcd=x;
	rw=0;
	rs=1;
	en=1;
	delay();
	en=0;
}

void lcd_init(){
	cmd_sr(0x38);
	delay();
	cmd_sr(0x01);
	delay();
	cmd_sr(0x0E);
	delay();
	cmd_sr(0x06);
	delay();
	cmd_sr(0x80);
	delay();
}
void timer_5s(){
	unsigned int i;
	for(i=0;i<76;i++){
		TF0=0;
		TH0=0x00;
		TL0=0x00;
		TR0=1;
		while(TF0==0);
		TR0=0;
		TF0=0;
	}
	return;
}
void entry() interrupt 0{
	unsigned int i;
	if(slots>0){
		cmd_sr(0x01);
		for(i=0; i<strlen(pwmsg); i++){
			data_sr(pwmsg[i]);
		}
	delay();

        if (pw1 == 1 && pw2 == 0 && pw3 == 0 && pw4 == 1) {
            slots--;

            doors = 0x01; delay();
            doors = 0x02; delay();
            doors = 0x04; delay();
            doors = 0x08; delay();

            timer_5s();

            doors = 0x04; delay();
            doors = 0x02; delay();
            doors = 0x01; delay();
        }
				else{
					cmd_sr(0x01);
		for(i=0; i<strlen(wrong); i++){
			data_sr(wrong[i]);
		}
    }
	}
}
void exit() interrupt 2{
	if(slots<5){
	slots++;
	doors = 0x08; delay();
        doors = 0x40; delay();
        doors = 0x20; delay();
        doors = 0x10; delay();

        timer_5s();

        doors = 0x20; delay();
        doors = 0x40; delay();
        doors = 0x80; delay();
    }
}
void main(){
	unsigned int i;
	pw1 = pw2 = pw3 = pw4 = 1;
	lcd_init();
	TMOD=0x01;
	EX0=1;
	IT0=1;
	EX1=1;
	IT1=1;
	EA=1;
	cmd_sr(0x01);
	delay();

	while(1){
		cmd_sr(0x01);
		delay();

		if(slots > 0){
			cmd_sr(0x80);  
			for(i = 0; i < strlen(str1); i++){
				data_sr(str1[i]);
			}

			cmd_sr(0xC0);
			for(i = 0; i < strlen(str2); i++){
				data_sr(str2[i]);
			}

			data_sr('0' + slots);  
		}
		else if(slots==0){
			cmd_sr(0x80);
			for(i = 0; i < strlen(msg); i++){
				data_sr(msg[i]);
			}
		}

		delay(); 
	}
}
